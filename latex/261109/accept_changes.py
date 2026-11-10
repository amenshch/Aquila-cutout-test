#!/usr/bin/env python3
"""
accept_changes.py -- accept all marked changes of a LaTeX file.

The text inside \\del{...} is removed together with the command, and the text inside \\add{...} is kept without the command.
Only the body of the document, after \\begin{document}, is processed, so the definitions of the two commands are kept.
Blank lines and spaces that a removal leaves doubled are merged into one.

Usage: python3 accept_changes.py <marked.tex> <clean.tex>
"""

import re
import sys


def closing(s, i):
    """The index of the brace that closes the one opened at s[i]."""
    depth = 0
    while i < len(s):
        c = s[i]
        if c == "\\":
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    sys.exit(" unbalanced braces")


def accept(body):
    out, i = [], 0
    while True:
        m = re.compile(r"\\(add|del)\{").search(body, i)
        if not m:
            out.append(body[i:])
            break
        out.append(body[i:m.start()])
        j = closing(body, m.end() - 1)
        i = j + 1
        if m.group(1) == "add":
            out.append(accept(body[m.end():j]))
        elif out and out[-1][-1:] in " \n":
            while i < len(body) and body[i] == " ":
                i += 1
    return "".join(out)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return
    s = open(sys.argv[1]).read()
    k = s.index(r"\begin{document}")
    body = accept(s[k:])
    body = re.sub(r"\n[ \t]*\n([ \t]*\n)+", "\n\n", body)
    body = re.sub(r"[ \t]+\n", "\n", body)
    open(sys.argv[2], "w").write(s[:k] + body)
    print(" written", sys.argv[2])


if __name__ == "__main__":
    main()
