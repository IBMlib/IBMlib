#!/usr/bin/env python
# -----------------------------------------------------------------
# replace all whitespaces (except "\n") with " "
# usage: cat test.bgm|regularize_blanks.py > test_only_spaces.bgm
# -----------------------------------------------------------------
import string
from   sys    import *
#
content = stdin.read()
for blank in string.whitespace:
    if blank == "\n":
        continue
    else:
        content = content.replace(blank, " ")
stdout.write(content)
#
