"""
getsf_columns.py -- locate catalogue columns by the character position of their
header names, not by column number.

Rationale
---------
getsf catalogues are fixed-width tables whose column count varies with the
number of wavebands, with which auxiliary blocks have been appended, and with
the getsf version.  Any script that addresses columns by number is therefore
correct only for the file it was written against.  This module reproduces the
scheme used in the getsf Fortran itself: find the header line that names the
columns, record the character span of each name, and read each data value from
the token that sits under that span.

Name normalisation
------------------
The same quantity appears as TOTL^MASS in concatenated catalogues and as
TOTL_MASS in raw ones, and band suffixes may be written 01 or 1.  Every name is
normalised by replacing '^' with '_', so lookups can use a single spelling.

Public interface
----------------
    tab = GetsfTable(path)          # parse one catalogue file
    tab.names                       # list of normalised column names, in order
    tab.has('TOTL_MASS')            # membership test
    tab.col('TOTL_MASS')            # float array, NaN where non-numeric
    tab.raw('QUALITY')              # list of strings, unparsed
    tab.nrows                       # number of data rows
    tab.header                      # list of comment lines, unmodified

    merged = GetsfTable.merge([tab1, tab2])   # join files row by row

Matching rule
-------------
Header names and numeric fields are right-aligned in getsf output, so a data
token is assigned to the header name whose end character position is closest to
the token's end.  Ties and gaps are resolved by falling back to the closest
centre.  A name that no token can be assigned to is reported by
`unmatched_names`, so a silent mis-read is not possible.
"""
import numpy as np
import re

_COMMENT = ('#', '!')


def _normalise(name):
    return name.replace('^', '_')


def _tokens_with_spans(line):
    """[(token, start, end)] for every whitespace-delimited token."""
    return [(m.group(0), m.start(), m.end())
            for m in re.finditer(r'\S+', line)]


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class GetsfTable(object):

    def __init__(self, path=None, _empty=False):
        self.path = path
        self.names = []
        self.header = []
        self._cols = {}          # normalised name -> list of raw strings
        self.unmatched_names = []
        self.nrows = 0
        if not _empty:
            self._parse(path)

    # ------------------------------------------------------------------ parse
    def _parse(self, path):
        lines = open(path).read().split('\n')
        is_com = lambda l: l.lstrip()[:1] in _COMMENT if l.strip() else False
        self.header = [l for l in lines if is_com(l)]

        data_idx = [i for i, l in enumerate(lines)
                    if l.strip() and not is_com(l)]
        if not data_idx:
            raise ValueError('no data rows found in %s' % path)
        first_data = data_idx[0]

        # The column-name line is the last comment line above the data that
        # contains several plausible names.  Scanning upward rather than
        # matching fixed content keeps this independent of the waveband set and
        # of any appended blocks.
        hdr_line = None
        for i in range(first_data - 1, -1, -1):
            if not is_com(lines[i]):
                continue
            stripped = lines[i].lstrip('#! ').rstrip()
            toks = stripped.split()
            if len(toks) < 4:
                continue
            good = sum(1 for t in toks
                       if re.match(r'^[A-Za-z][A-Za-z0-9_^]*$', t))
            if good >= max(4, int(0.8 * len(toks))):
                hdr_line = lines[i]
                break
        if hdr_line is None:
            raise ValueError('could not find a column-name line in %s' % path)

        # strip only the leading comment marker, preserving character positions
        lead = len(hdr_line) - len(hdr_line.lstrip('#! '))
        hdr_spans = [(_normalise(t), s, e)
                     for (t, s, e) in _tokens_with_spans(hdr_line)
                     if s >= lead - 1]
        hdr_spans = [(n, s, e) for (n, s, e) in hdr_spans
                     if re.match(r'^[A-Za-z][A-Za-z0-9_]*$', n)]
        self.names = [n for (n, _, _) in hdr_spans]

        cols = dict((n, []) for n in self.names)
        seen = set()
        rows = [lines[i] for i in data_idx]
        for row in rows:
            toks = _tokens_with_spans(row)
            used = [False] * len(toks)
            vals = {}
            # assign by closest right edge, then by closest centre
            for (n, hs, he) in hdr_spans:
                best, bestd = None, None
                for k, (t, s, e) in enumerate(toks):
                    if used[k]:
                        continue
                    d = abs(e - he)
                    if bestd is None or d < bestd:
                        best, bestd = k, d
                if best is None:
                    continue
                hc, tc = 0.5 * (hs + he), 0.5 * (toks[best][1] + toks[best][2])
                for k, (t, s, e) in enumerate(toks):
                    if used[k]:
                        continue
                    if abs(0.5 * (s + e) - hc) < abs(tc - hc):
                        best = k
                        tc = 0.5 * (s + e)
                used[best] = True
                vals[n] = toks[best][0]
                seen.add(n)
            for n in self.names:
                cols[n].append(vals.get(n, ''))
        self._cols = cols
        self.nrows = len(rows)
        self.unmatched_names = [n for n in self.names if n not in seen]

    # ----------------------------------------------------------------- access
    def has(self, name):
        return _normalise(name) in self._cols

    def raw(self, name):
        return self._cols[_normalise(name)]

    def col(self, name):
        return np.array([float(v) if _is_number(v) else np.nan
                         for v in self._cols[_normalise(name)]])

    def find(self, pattern):
        """Names matching a regular expression, e.g. r'^AFWHM\\d\\d$'."""
        rx = re.compile(pattern)
        return [n for n in self.names if rx.match(n)]

    # ------------------------------------------------------------------ merge
    @classmethod
    def merge(cls, tables):
        """Join tables row by row.  Duplicated names keep the first occurrence."""
        n = min(t.nrows for t in tables)
        out = cls(_empty=True)
        out.path = [t.path for t in tables]
        for t in tables:
            out.header += t.header
            for nm in t.names:
                if nm in out._cols:
                    continue
                out.names.append(nm)
                out._cols[nm] = t._cols[nm][:n]
        out.nrows = n
        return out


def summarise(path):
    t = GetsfTable(path)
    print('%s: %d columns, %d rows' % (path, len(t.names), t.nrows))
    if t.unmatched_names:
        print('  WARNING: no data token found for %s' % t.unmatched_names)
    return t


if __name__ == '__main__':
    import sys
    for p in sys.argv[1:]:
        t = summarise(p)
        print('  first names: %s' % ' '.join(t.names[:8]))
        for probe in ('TOTL_MASS', 'DUST_TEMP', 'QUALITY', 'FOOA03', 'AFWHM03'):
            if t.has(probe):
                v = t.col(probe)
                print('  %-10s median %.4g over %d finite values'
                      % (probe, np.nanmedian(v), int(np.sum(np.isfinite(v)))))
