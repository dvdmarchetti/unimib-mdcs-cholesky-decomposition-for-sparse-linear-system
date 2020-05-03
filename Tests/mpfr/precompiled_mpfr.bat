@ECHO OFF
IF NOT "%~1"=="" GOTO lbl_ok
TYPE _internal\precompiled-help
EXIT /B 1
:lbl_ok
COPY _internal\"%~1"-mpn.lib mpfr\mpn.lib > NUL
COPY _internal\"%~1"-gmp.lib mpfr\gmp.lib > NUL
COPY _internal\"%~1"-mpfr.lib mpfr\mpfr.lib > NUL
