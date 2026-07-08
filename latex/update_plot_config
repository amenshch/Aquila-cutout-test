#!/usr/bin/env bash
# update_plot_cfgs — Regenerate +densities.cfg, +surfdens.cfg, +tempers.cfg
# in the current tBE_j directory to reflect all existing zk subdirectories.
# Run from inside a tBE_j directory.

installpath=${GETSF_HOME/_REVISIONS/}
source $installpath'/comfun'
verb='-verb2'; logfile='/dev/null'; cwavenum=''
set -euo pipefail
unset -f find 2>/dev/null || true

# ---------------------------------------------------------------------------
# Plot axis ranges (edit as needed)
# ---------------------------------------------------------------------------
DENS_XMIN='10.0'
DENS_XMAX='1.0e7'
DENS_YMIN='1.0e-23'
DENS_YMAX='1.0e-17'

SURF_XMIN='10.0'
SURF_XMAX='1.0e7'
SURF_YMIN='1.0e+20'
SURF_YMAX='2.5e+23'

TEMP_XMIN='10.0'
TEMP_XMAX='1.0e7'
TEMP_YMIN='6.0'
TEMP_YMAX='15.0'

# ---------------------------------------------------------------------------
# Discover existing zk subdirs (numeric, 2-digit), sorted
# ---------------------------------------------------------------------------
mapfile -t ZK_DIRS < <(\find . -maxdepth 1 -type d -name '[0-9][0-9]' | \
                        sed 's|^\./||' | sort)

if [[ ${#ZK_DIRS[@]} -eq 0 ]]; then
    echo "No numeric subdirectories found in $(pwd). Aborting." >&2
    exit 1
fi

echo; echo "Found zk dirs: ${ZK_DIRS[*]}"

# ---------------------------------------------------------------------------
# Blank the var field in a block (for all dirs after the first).
# Matches the var line by key and replaces the leading quoted string with ''.
# ---------------------------------------------------------------------------
blank_var() {
    # Input: multiline block on stdin
    # Only the line containing "var    :" is modified
    awk '{
        if (/var[[:space:]]+:/) {
            # Replace leading quoted token (possibly containing |) with empty quotes
            sub(/^[[:space:]]*'"'"'[^'"'"']*'"'"'/, "'"'"''"'"'")
        }
        print
    }'
}

# ---------------------------------------------------------------------------
# Generic cfg writer
# $1: output filename
# $2: header string
# $3: footer string
# $4+: block template files (DIRNUM substituted per directory)
# ---------------------------------------------------------------------------
write_cfg() {
    local outfile="$1"
    local header="$2"
    local footer="$3"
    shift 3
    local -a block_files=("$@")

    {
        printf '%s\n' "$header"

        local first_dir=1
        for zk in "${ZK_DIRS[@]}"; do
            for btpl in "${block_files[@]}"; do
                local block
                block=$(sed "s|DIRNUM|${zk}|g" "$btpl")

                if [[ $first_dir -eq 0 ]]; then
                    block=$(printf '%s\n' "$block" | blank_var)
                fi

                printf '%s\n' "$block"
                printf ';\n'
            done
            first_dir=0
        done

        printf '%s\n' "$footer"
    } > "$outfile"

    echo "  Written: $outfile"
}

# ---------------------------------------------------------------------------
# Write block templates to temp files
# ---------------------------------------------------------------------------
TMPDIR_LOCAL=$(mktemp -d)
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

# --- +densities.cfg blocks ---

cat > "$TMPDIR_LOCAL/dens_blk1.tpl" << 'TPLEOF'
'DIRNUM/mc.denstemps.tab' ............................ fname  : data file name
'RADIUS(au)' ..................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'RELZONESIZE' .................................... nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 1.00001 ......................................... ascale : {1.0} argument scaling factor
 1.0e-20 ......................................... fscale : {1.0} function scaling factor
'|F5|d|F10|r|F5|/|F10|r' ......................... var    : variable name (none if blank)
'gra' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.0 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

cat > "$TMPDIR_LOCAL/dens_blk2.tpl" << 'TPLEOF'
'DIRNUM/mc.denstemps.tab' ............................ fname  : data file name
'RADIUS(au)' ..................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'DENS(g/cm^3)' ................................... nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 1.00001 ......................................... ascale : {1.0} argument scaling factor
 1.00002 ......................................... fscale : {1.0} function scaling factor
'|F5|r' .......................................... var    : variable name (none if blank)
'red' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.5 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

# --- +surfdens.cfg blocks ---

cat > "$TMPDIR_LOCAL/surf_blk1.tpl" << 'TPLEOF'
'DIRNUM/mc.surfdens.pro.bg.tab' ...................... fname  : data file name
'RADIUS(au)' ..................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'PROFILEBG' ...................................... nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 1.00001 ......................................... ascale : {1.0} argument scaling factor
 1.00002 ......................................... fscale : {1.0} function scaling factor
'|F5|s|BF0|bgl' .................................. var    : variable name (none if blank)
'blu' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.5 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

cat > "$TMPDIR_LOCAL/surf_blk2.tpl" << 'TPLEOF'
'DIRNUM/mc.surfdens.pro.tab' ......................... fname  : data file name
'RADIUS(au)' ..................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'PROFILE' ........................................ nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 1.00001 ......................................... ascale : {1.0} argument scaling factor
 1.00002 ......................................... fscale : {1.0} function scaling factor
'|F5|s' .......................................... var    : variable name (none if blank)
'red' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.5 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

cat > "$TMPDIR_LOCAL/surf_blk3.tpl" << 'TPLEOF'
'DIRNUM/mc.surfdens.pro.tab' ......................... fname  : data file name
'RADIUS(au)' ..................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'PROFILEBS' ...................................... nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 1.00001 ......................................... ascale : {1.0} argument scaling factor
 1.00002 ......................................... fscale : {1.0} function scaling factor
'|F5|s|BF0|bs' ................................... var    : variable name (none if blank)
'gre' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.5 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

cat > "$TMPDIR_LOCAL/surf_blk4.tpl" << 'TPLEOF'
'DIRNUM/mc.surfdens.bsl.pro.dat' ..................... fname  : data file name
'ARCSEC' ......................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'VALUE' .......................................... nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 260 ............................................. ascale : {1.0} argument scaling factor
 1.00002 ......................................... fscale : {1.0} function scaling factor
'|F5|s|BF0|bsl' .................................. var    : variable name (none if blank)
'mag' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.5 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

# --- +tempers.cfg blocks ---

cat > "$TMPDIR_LOCAL/temp_blk1.tpl" << 'TPLEOF'
'DIRNUM/mc.denstemps.tab' ............................ fname  : data file name
'RADIUS(au)' ..................................... nxc    : column name or number for argument X; if a 2nd column name follows '%', then X1/X2 is plotted
'TEMPERATURE' .................................... nyc    : column name or number for function Y; if a 2nd column name follows '%', then Y1/Y2 is plotted
 0 ............................................... nexpc  : column name or number for errX+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... nexmc  : column name or number for errX- (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neypc  : column name or number for errY+ (< 0: as percentage, ?? or 0: no errors)
 0 ............................................... neymc  : column name or number for errY- (< 0: as percentage, ?? or 0: no errors)
 1.00001 ......................................... ascale : {1.0} argument scaling factor
 1.00002 ......................................... fscale : {1.0} function scaling factor
'|F10|T|BF0|d' ................................... var    : variable name (none if blank)
'red' ............................................ color  : {red gre blu cya mag bro ora gol dre dgr dbl dcy dma dbr dor dgo bla gra#; auto}
 0 ............................................... ndsh   : {2-4-8-16-} dash styles (0: solid lines)   {xred xgre xblu xora xpur xcya xmag xtea xyel xpin xlim}
 1 ............................................... lwid   : {-8-6-4-2-1} line width (0: no lines)      {xlav xbro xbei xmar xmin xoli xapr xnav xgra xwhi xbla}
 0.5 ............................................. ssze   : {-0.5-1.5-} symbol size (0.0: no symbols)
'circle' ......................................... symbol : {none circle square triangle trianglei diamond minus}
TPLEOF

# ---------------------------------------------------------------------------
# Headers and footers
# ---------------------------------------------------------------------------

DENS_HEADER=';
; Parameters for 1D plotting in PostScript (plot):        ital |F10| greek |F5| deflt |F0| sun |BF18|* deg |F18|o asec |F18|"
;
'"'"'Total Density Distributions'"'"' .................... phdr   : plot title (none if blank or starts with ??)
'"'"'Radial distance |F10|r|F0| (AU)'"'"' ................ xlab   : label of the main horizontal axis (none if blank or starts with ??)
'"'"'Density |F5|r|F0|(|F10|r|F0|) (g|B| |N|cm|S|-3|N|)'"'"' ......... ylab   : label of the main vertical axis (none if blank or starts with ??)
'"'"'Angular distance |F5|q|F0| (|F18|"|N|)'"'"' ......... xlabu  : label of the upper horizontal axis (none if blank or starts with ??)
'"'"' '"'"' '"'"'Density |F5|r|F0| (g|B| |N|cm|S|-3|N|)'"'"' ..... ylabr  : label of the right vertical axis (none if blank or starts with ??)
'"'"'portrait'"'"' ....................................... orient : {portrait landscape} page orientation
10.00001 ......................................... xsize  : horizontal size of the plot (~12 or ~24 cm)
10.00002 ......................................... ysize  : vertical size of the plot (~12 cm)
 4.5 ............................................. plt1x0 : x-coordinate of the lower left corner (~3.3 or ~2.3 cm)
17.6 ............................................. plt1y0 : y-coordinate of the lower left corner (~14.5 or ~5.5 cm)
'"'"'loglog'"'"' ......................................... axtype : {linlin linlog loglin loglog} axes types
'"'"'nono'"'"' ........................................... fnorm  : {norm nono} normalize function to 1.0 or not
 '"${DENS_XMIN}"' ............................................. xmin   : minimum value of argument X to plot (0.0: automatic)
 '"${DENS_XMAX}"' ........................................... xmax   : maximum value of argument X to plot (0.0: automatic)
 0.003846154 1.0 0.0 ............................. upper  : coeffs a b c in a*X^b+c for upper axis (0.0 * *: none)
 '"${DENS_YMIN}"' ......................................... ymin   : minimum value of function Y to plot (0.0: automatic)
 '"${DENS_YMAX}"' ......................................... ymax   : maximum value of function Y to plot (0.0: automatic)
 0.0 1.0 0.0 ..................................... right  : coeffs d e f in d*Y^e+f for right axis (0.0 * *: none)
 2 ............................................... allxl  : {-1,0,1,2} plot all labels/annotations on the X side
 2 ............................................... allyl  : {0,1,2} plot all labels/annotations on the Y side
'"'"'allxh'"'"' .......................................... xhlp   : {allxh major none} help tick lines normal to X axis
'"'"'allyh'"'"' .......................................... yhlp   : {allyh major none} help tick lines normal to Y axis
'"'"'manual'"'"' ......................................... colors : {auto manual} choose colors automatically or not
0 ................................................ nfun   : number of functions to plot (0: plot all)
;
; Parameters for each variable (up to 999 functions):          sub |B| sup |S| norm |N| times |F18|x .ge. |F18|>= prop |F18|~
;'

DENS_FOOTER='; Text lines to be overlaid on plot (as many lines as needed):
;
'"'"''"'"' 3.0 3.0 1.0 ................................... textxy : {annotation X Y sizefactor} (cm) to overlay on plot'

SURF_HEADER=';
; Parameters for 1D plotting in PostScript (plot):        ital |F10| greek |F5| deflt |F0| sun |BF18|* deg |F18|o asec |F18|"
;
'"'"''"'"' ............................................... phdr   : plot title (none if blank or starts with ??)
'"'"'Radial distance |F10|r|F0| (AU)'"'"' ................ xlab   : label of the main horizontal axis (none if blank or starts with ??)
'"'"'Surface density |F5|s|BF0|H|B|2|NF0|(|F10|r|F0|) (cm|S|-2|N|)'"'"' ......... ylab   : label of the main vertical axis (none if blank or starts with ??)
'"'"'Angular distance |F5|q|F0| (|F18|"|N|)'"'"' ......... xlabu  : label of the upper horizontal axis (none if blank or starts with ??)
'"'"' '"'"' .............................................. ylabr  : label of the right vertical axis (none if blank or starts with ??)
'"'"'portrait'"'"' ....................................... orient : {portrait landscape} page orientation
10.00001 ......................................... xsize  : horizontal size of the plot (~12 or ~24 cm)
10.00002 ......................................... ysize  : vertical size of the plot (~12 cm)
 4.5 ............................................. plt1x0 : x-coordinate of the lower left corner (~3.3 or ~2.3 cm)
17.6 ............................................. plt1y0 : y-coordinate of the lower left corner (~14.5 or ~5.5 cm)
'"'"'loglog'"'"' ......................................... axtype : {linlin linlog loglin loglog} axes types
'"'"'nono'"'"' ........................................... fnorm  : {norm nono} normalize function to 1.0 or not
 '"${SURF_XMIN}"' ............................................. xmin   : minimum value of argument X to plot (0.0: automatic)
 '"${SURF_XMAX}"' ........................................... xmax   : maximum value of argument X to plot (0.0: automatic)
 0.003846154 1.0 0.0 ............................. upper  : coeffs a b c in a*X^b+c for upper axis (0.0 * *: none)
 '"${SURF_YMIN}"' ......................................... ymin   : minimum value of function Y to plot (0.0: automatic)
 '"${SURF_YMAX}"' ......................................... ymax   : maximum value of function Y to plot (0.0: automatic)
 0.0 1.0 0.0 ..................................... right  : coeffs d e f in d*Y^e+f for right axis (0.0 * *: none)
 2 ............................................... allxl  : {-1,0,1,2} plot all labels/annotations on the X side
 2 ............................................... allyl  : {0,1,2} plot all labels/annotations on the Y side
'"'"'allxh'"'"' .......................................... xhlp   : {allxh major none} help tick lines normal to X axis
'"'"'allyh'"'"' .......................................... yhlp   : {allyh major none} help tick lines normal to Y axis
'"'"'manual'"'"' ......................................... colors : {auto manual} choose colors automatically or not
0 ................................................ nfun   : number of functions to plot (0: plot all)
;
; Parameters for each variable (up to 999 functions):          sub |B| sup |S| norm |N| times |F18|x .ge. |F18|>= prop |F18|~
;'

SURF_FOOTER='; Text lines to be overlaid on plot (as many lines as needed):
;
'"'"''"'"' 3.0 3.0 1.0 ................................... textxy : {annotation X Y sizefactor} (cm) to overlay on plot'

TEMP_HEADER=';
; Parameters for 1D plotting in PostScript (plot):        ital |F10| greek |F5| deflt |F0| sun |BF18|* deg |F18|o asec |F18|"
;
'"'"'Temperature Distributions'"'"' ...................... phdr   : plot title (none if blank or starts with ??)
'"'"'Radial distance |F10|r|F0| (AU)'"'"' ................ xlab   : label of the main horizontal axis (none if blank or starts with ??)
'"'"'Dust temperature |F10|T|BF0|d|NF0|(|F10|r|F0|) (K)'"'"' ........... ylab   : label of the main vertical axis (none if blank or starts with ??)
'"'"'Angular distance |F5|q|F0| (|F18|"|N|)'"'"' ......... xlabu  : label of the upper horizontal axis (none if blank or starts with ??)
'"'"' '"'"' '"'"'Density |F5|r|F0| (g|B| |N|cm|S|-3|N|)'"'"' ..... ylabr  : label of the right vertical axis (none if blank or starts with ??)
'"'"'portrait'"'"' ....................................... orient : {portrait landscape} page orientation
10.00001 ......................................... xsize  : horizontal size of the plot (~12 or ~24 cm)
10.00002 ......................................... ysize  : vertical size of the plot (~12 cm)
 4.5 ............................................. plt1x0 : x-coordinate of the lower left corner (~3.3 or ~2.3 cm)
17.6 ............................................. plt1y0 : y-coordinate of the lower left corner (~14.5 or ~5.5 cm)
'"'"'loglog'"'"' ......................................... axtype : {linlin linlog loglin loglog} axes types
'"'"'nono'"'"' ........................................... fnorm  : {norm nono} normalize function to 1.0 or not
 '"${TEMP_XMIN}"' ............................................. xmin   : minimum value of argument X to plot (0.0: automatic)
 '"${TEMP_XMAX}"' ........................................... xmax   : maximum value of argument X to plot (0.0: automatic)
 0.003846154 1.0 0.0 ............................. upper  : coeffs a b c in a*X^b+c for upper axis (0.0 * *: none)
 '"${TEMP_YMIN}"' ......................................... ymin   : minimum value of function Y to plot (0.0: automatic)
 '"${TEMP_YMAX}"' ......................................... ymax   : maximum value of function Y to plot (0.0: automatic)
 0.0 1.0 0.0 ..................................... right  : coeffs d e f in d*Y^e+f for right axis (0.0 * *: none)
 2 ............................................... allxl  : {-1,0,1,2} plot all labels/annotations on the X side
 2 ............................................... allyl  : {0,1,2} plot all labels/annotations on the Y side
'"'"'allxh'"'"' .......................................... xhlp   : {allxh major none} help tick lines normal to X axis
'"'"'allyh'"'"' .......................................... yhlp   : {allyh major none} help tick lines normal to Y axis
'"'"'manual'"'"' ......................................... colors : {auto manual} choose colors automatically or not
0 ................................................ nfun   : number of functions to plot (0: plot all)
;
; Parameters for each variable (up to 999 functions):          sub |B| sup |S| norm |N| times |F18|x .ge. |F18|>= prop |F18|~
;'

TEMP_FOOTER='; Text lines to be overlaid on plot (as many lines as needed):
;
'"'"''"'"' 3.0 3.0 1.0 ................................... textxy : {annotation X Y sizefactor} (cm) to overlay on plot'

# ---------------------------------------------------------------------------
# Generate all three cfg files
# ---------------------------------------------------------------------------

echo "Updating cfg files in $(pwd) ..."
write_cfg '+densities.cfg' "$DENS_HEADER" "$DENS_FOOTER" \
    "$TMPDIR_LOCAL/dens_blk1.tpl" "$TMPDIR_LOCAL/dens_blk2.tpl"
write_cfg '+surfdens.cfg'  "$SURF_HEADER" "$SURF_FOOTER" \
    "$TMPDIR_LOCAL/surf_blk1.tpl" "$TMPDIR_LOCAL/surf_blk2.tpl" \
    "$TMPDIR_LOCAL/surf_blk3.tpl" "$TMPDIR_LOCAL/surf_blk4.tpl"
write_cfg '+tempers.cfg'   "$TEMP_HEADER" "$TEMP_FOOTER" \
    "$TMPDIR_LOCAL/temp_blk1.tpl"

('plot' '+densities' $verb; echo $? > '.rc') | 'tehi'
('plot' '+surfdens' $verb; echo $? > '.rc') | 'tehi'
('plot' '+tempers' $verb; echo $? > '.rc') | 'tehi'

echo "Done."
