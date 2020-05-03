@echo off
IF NOT "%~1"=="" GOTO lbl_1_end
type _internal\build-bat-help
EXIT /B 1
:lbl_1_end
IF "%~3"=="" GOTO lbl_2_end
echo Too many parameters specified
echo Did you enclose parameters to be passed to the compiler in double quotes?
EXIT /B 1
:lbl_2_end
CHDIR out
IF EXIST * DEL /F /Q *
CHDIR ..\_tmp
IF EXIST * DEL /F /Q *
CHDIR ..
COPY src\* _tmp > NUL 2> NUL
IF "%~1"=="msvc" GOTO lbl_3_msvc
GOTO lbl_3___error
:lbl_3_msvc
CHDIR _tmp
cl /nologo /EHsc /I. /c %~2 *.cpp  >> ../log.txt 2>&1
IF NOT ERRORLEVEL 1 GOTO lbl_4
echo Error while compiling (see ../log.txt for more info)
CHDIR ..
EXIT /B 1
:lbl_4
lib /out:libmpalglib.lib *.obj >> ../log.txt 2>> ../log.txt 2>&1
IF NOT ERRORLEVEL 1 GOTO lbl_5
echo Error while running LIB (see ../log.txt for more info)
CHDIR ..
EXIT /B 1
:lbl_5
CHDIR ..
COPY _tmp\libmpalglib.lib out > NUL 2> NUL
COPY src\*.h out > NUL 2> NUL
COPY mpfr\mpn.lib out\mpn.lib > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_6
echo Error copying mpn static library. Did you ran precompiled_mpfr?
EXIT /B 1
:lbl_6
COPY mpfr\gmp.lib out\gmp.lib > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_7
echo Error copying gmp static library. Did you ran precompiled_mpfr?
EXIT /B 1
:lbl_7
COPY mpfr\mpfr.lib out\mpfr.lib > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_8
echo Error copying mpfr static library. Did you ran precompiled_mpfr?
EXIT /B 1
:lbl_8
CHDIR _tmp
IF EXIST * DEL /F /Q *
CHDIR ..
GOTO lbl_3___end
:lbl_3___error
echo unknown compiler
EXIT /B 1
:lbl_3___end
