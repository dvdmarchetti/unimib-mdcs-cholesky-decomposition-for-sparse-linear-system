@echo off
IF NOT "%~1"=="" GOTO lbl_1_end
type _internal\check-bat-help
EXIT /B 1
:lbl_1_end
IF "%~4"=="" GOTO lbl_2_end
echo Too many parameters specified
echo Did you enclose parameters to be passed to the compiler in double quotes?
EXIT /B 1
:lbl_2_end
IF NOT "%~2"=="all" GOTO lbl_3_end
pushd .
CALL .\check  "%~1"  "ablas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_4
EXIT /B 1
:lbl_4
pushd .
CALL .\check  "%~1"  "ablasf_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_5
EXIT /B 1
:lbl_5
pushd .
CALL .\check  "%~1"  "ortfac_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_6
EXIT /B 1
:lbl_6
pushd .
CALL .\check  "%~1"  "hblas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_7
EXIT /B 1
:lbl_7
pushd .
CALL .\check  "%~1"  "reflections_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_8
EXIT /B 1
:lbl_8
pushd .
CALL .\check  "%~1"  "creflections_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_9
EXIT /B 1
:lbl_9
pushd .
CALL .\check  "%~1"  "sblas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_10
EXIT /B 1
:lbl_10
pushd .
CALL .\check  "%~1"  "evd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_11
EXIT /B 1
:lbl_11
pushd .
CALL .\check  "%~1"  "blas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_12
EXIT /B 1
:lbl_12
pushd .
CALL .\check  "%~1"  "rotations_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_13
EXIT /B 1
:lbl_13
pushd .
CALL .\check  "%~1"  "hsschur_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_14
EXIT /B 1
:lbl_14
pushd .
CALL .\check  "%~1"  "matgen_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_15
EXIT /B 1
:lbl_15
pushd .
CALL .\check  "%~1"  "hqrnd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_16
EXIT /B 1
:lbl_16
pushd .
CALL .\check  "%~1"  "trfac_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_17
EXIT /B 1
:lbl_17
pushd .
CALL .\check  "%~1"  "rcond_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_18
EXIT /B 1
:lbl_18
pushd .
CALL .\check  "%~1"  "trlinsolve_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_19
EXIT /B 1
:lbl_19
pushd .
CALL .\check  "%~1"  "safesolve_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_20
EXIT /B 1
:lbl_20
pushd .
CALL .\check  "%~1"  "matinv_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_21
EXIT /B 1
:lbl_21
pushd .
CALL .\check  "%~1"  "bdsvd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_22
EXIT /B 1
:lbl_22
pushd .
CALL .\check  "%~1"  "svd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_23
EXIT /B 1
:lbl_23
pushd .
CALL .\check  "%~1"  "densesolver_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_24
EXIT /B 1
:lbl_24
pushd .
CALL .\check  "%~1"  "xblas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_25
EXIT /B 1
:lbl_25
pushd .
CALL .\check  "%~1"  "minlbfgs_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_26
EXIT /B 1
:lbl_26
pushd .
CALL .\check  "%~1"  "linmin_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_27
EXIT /B 1
:lbl_27
pushd .
CALL .\check  "%~1"  "odesolver_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_28
EXIT /B 1
:lbl_28
pushd .
CALL .\check  "%~1"  "conv_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_29
EXIT /B 1
:lbl_29
pushd .
CALL .\check  "%~1"  "ftbase_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_30
EXIT /B 1
:lbl_30
pushd .
CALL .\check  "%~1"  "fft_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_31
EXIT /B 1
:lbl_31
pushd .
CALL .\check  "%~1"  "corr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_32
EXIT /B 1
:lbl_32
pushd .
CALL .\check  "%~1"  "fht_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_33
EXIT /B 1
:lbl_33
pushd .
CALL .\check  "%~1"  "autogk_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_34
EXIT /B 1
:lbl_34
pushd .
CALL .\check  "%~1"  "tsort_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_35
EXIT /B 1
:lbl_35
pushd .
CALL .\check  "%~1"  "gammafunc_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_36
EXIT /B 1
:lbl_36
pushd .
CALL .\check  "%~1"  "gq_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_37
EXIT /B 1
:lbl_37
pushd .
CALL .\check  "%~1"  "gkq_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_38
EXIT /B 1
:lbl_38
pushd .
CALL .\check  "%~1"  "lsfit_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_39
EXIT /B 1
:lbl_39
pushd .
CALL .\check  "%~1"  "minlm_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_40
EXIT /B 1
:lbl_40
pushd .
CALL .\check  "%~1"  "polint_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_41
EXIT /B 1
:lbl_41
pushd .
CALL .\check  "%~1"  "ratinterpolation_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_42
EXIT /B 1
:lbl_42
pushd .
CALL .\check  "%~1"  "ratint_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_43
EXIT /B 1
:lbl_43
pushd .
CALL .\check  "%~1"  "apserv_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_44
EXIT /B 1
:lbl_44
pushd .
CALL .\check  "%~1"  "spline2d_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_45
EXIT /B 1
:lbl_45
pushd .
CALL .\check  "%~1"  "spline3_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_46
EXIT /B 1
:lbl_46
pushd .
CALL .\check  "%~1"  "spline1d_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_47
EXIT /B 1
:lbl_47
pushd .
CALL .\check  "%~1"  "idwint_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_48
EXIT /B 1
:lbl_48
pushd .
CALL .\check  "%~1"  "nearestneighbor_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_49
EXIT /B 1
:lbl_49
pushd .
CALL .\check  "%~1"  "pspline_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_50
EXIT /B 1
:lbl_50
pushd .
CALL .\check  "%~1"  "matdet_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_51
EXIT /B 1
:lbl_51
pushd .
CALL .\check  "%~1"  "sdet_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_52
EXIT /B 1
:lbl_52
pushd .
CALL .\check  "%~1"  "ldlt_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_53
EXIT /B 1
:lbl_53
pushd .
CALL .\check  "%~1"  "spdgevd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_54
EXIT /B 1
:lbl_54
pushd .
CALL .\check  "%~1"  "sinverse_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_55
EXIT /B 1
:lbl_55
pushd .
CALL .\check  "%~1"  "inverseupdate_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_56
EXIT /B 1
:lbl_56
pushd .
CALL .\check  "%~1"  "srcond_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_57
EXIT /B 1
:lbl_57
pushd .
CALL .\check  "%~1"  "ssolve_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_58
EXIT /B 1
:lbl_58
pushd .
CALL .\check  "%~1"  "estnorm_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_59
EXIT /B 1
:lbl_59
pushd .
CALL .\check  "%~1"  "schur_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_60
EXIT /B 1
:lbl_60
pushd .
CALL .\check  "%~1"  "mincg_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_61
EXIT /B 1
:lbl_61
pushd .
CALL .\check  "%~1"  "minasa_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_62
EXIT /B 1
:lbl_62
EXIT /B 0
:lbl_3_end
IF NOT "%~2"=="all_silent" GOTO lbl_63_end
pushd .
CALL .\check  "%~1"  "ablas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_64
EXIT /B 1
:lbl_64
pushd .
CALL .\check  "%~1"  "ablasf_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_65
EXIT /B 1
:lbl_65
pushd .
CALL .\check  "%~1"  "ortfac_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_66
EXIT /B 1
:lbl_66
pushd .
CALL .\check  "%~1"  "hblas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_67
EXIT /B 1
:lbl_67
pushd .
CALL .\check  "%~1"  "reflections_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_68
EXIT /B 1
:lbl_68
pushd .
CALL .\check  "%~1"  "creflections_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_69
EXIT /B 1
:lbl_69
pushd .
CALL .\check  "%~1"  "sblas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_70
EXIT /B 1
:lbl_70
pushd .
CALL .\check  "%~1"  "evd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_71
EXIT /B 1
:lbl_71
pushd .
CALL .\check  "%~1"  "blas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_72
EXIT /B 1
:lbl_72
pushd .
CALL .\check  "%~1"  "rotations_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_73
EXIT /B 1
:lbl_73
pushd .
CALL .\check  "%~1"  "hsschur_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_74
EXIT /B 1
:lbl_74
pushd .
CALL .\check  "%~1"  "matgen_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_75
EXIT /B 1
:lbl_75
pushd .
CALL .\check  "%~1"  "hqrnd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_76
EXIT /B 1
:lbl_76
pushd .
CALL .\check  "%~1"  "trfac_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_77
EXIT /B 1
:lbl_77
pushd .
CALL .\check  "%~1"  "rcond_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_78
EXIT /B 1
:lbl_78
pushd .
CALL .\check  "%~1"  "trlinsolve_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_79
EXIT /B 1
:lbl_79
pushd .
CALL .\check  "%~1"  "safesolve_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_80
EXIT /B 1
:lbl_80
pushd .
CALL .\check  "%~1"  "matinv_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_81
EXIT /B 1
:lbl_81
pushd .
CALL .\check  "%~1"  "bdsvd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_82
EXIT /B 1
:lbl_82
pushd .
CALL .\check  "%~1"  "svd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_83
EXIT /B 1
:lbl_83
pushd .
CALL .\check  "%~1"  "densesolver_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_84
EXIT /B 1
:lbl_84
pushd .
CALL .\check  "%~1"  "xblas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_85
EXIT /B 1
:lbl_85
pushd .
CALL .\check  "%~1"  "minlbfgs_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_86
EXIT /B 1
:lbl_86
pushd .
CALL .\check  "%~1"  "linmin_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_87
EXIT /B 1
:lbl_87
pushd .
CALL .\check  "%~1"  "odesolver_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_88
EXIT /B 1
:lbl_88
pushd .
CALL .\check  "%~1"  "conv_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_89
EXIT /B 1
:lbl_89
pushd .
CALL .\check  "%~1"  "ftbase_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_90
EXIT /B 1
:lbl_90
pushd .
CALL .\check  "%~1"  "fft_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_91
EXIT /B 1
:lbl_91
pushd .
CALL .\check  "%~1"  "corr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_92
EXIT /B 1
:lbl_92
pushd .
CALL .\check  "%~1"  "fht_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_93
EXIT /B 1
:lbl_93
pushd .
CALL .\check  "%~1"  "autogk_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_94
EXIT /B 1
:lbl_94
pushd .
CALL .\check  "%~1"  "tsort_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_95
EXIT /B 1
:lbl_95
pushd .
CALL .\check  "%~1"  "gammafunc_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_96
EXIT /B 1
:lbl_96
pushd .
CALL .\check  "%~1"  "gq_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_97
EXIT /B 1
:lbl_97
pushd .
CALL .\check  "%~1"  "gkq_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_98
EXIT /B 1
:lbl_98
pushd .
CALL .\check  "%~1"  "lsfit_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_99
EXIT /B 1
:lbl_99
pushd .
CALL .\check  "%~1"  "minlm_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_100
EXIT /B 1
:lbl_100
pushd .
CALL .\check  "%~1"  "polint_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_101
EXIT /B 1
:lbl_101
pushd .
CALL .\check  "%~1"  "ratinterpolation_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_102
EXIT /B 1
:lbl_102
pushd .
CALL .\check  "%~1"  "ratint_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_103
EXIT /B 1
:lbl_103
pushd .
CALL .\check  "%~1"  "apserv_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_104
EXIT /B 1
:lbl_104
pushd .
CALL .\check  "%~1"  "spline2d_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_105
EXIT /B 1
:lbl_105
pushd .
CALL .\check  "%~1"  "spline3_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_106
EXIT /B 1
:lbl_106
pushd .
CALL .\check  "%~1"  "spline1d_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_107
EXIT /B 1
:lbl_107
pushd .
CALL .\check  "%~1"  "idwint_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_108
EXIT /B 1
:lbl_108
pushd .
CALL .\check  "%~1"  "nearestneighbor_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_109
EXIT /B 1
:lbl_109
pushd .
CALL .\check  "%~1"  "pspline_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_110
EXIT /B 1
:lbl_110
pushd .
CALL .\check  "%~1"  "matdet_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_111
EXIT /B 1
:lbl_111
pushd .
CALL .\check  "%~1"  "sdet_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_112
EXIT /B 1
:lbl_112
pushd .
CALL .\check  "%~1"  "ldlt_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_113
EXIT /B 1
:lbl_113
pushd .
CALL .\check  "%~1"  "spdgevd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_114
EXIT /B 1
:lbl_114
pushd .
CALL .\check  "%~1"  "sinverse_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_115
EXIT /B 1
:lbl_115
pushd .
CALL .\check  "%~1"  "inverseupdate_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_116
EXIT /B 1
:lbl_116
pushd .
CALL .\check  "%~1"  "srcond_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_117
EXIT /B 1
:lbl_117
pushd .
CALL .\check  "%~1"  "ssolve_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_118
EXIT /B 1
:lbl_118
pushd .
CALL .\check  "%~1"  "estnorm_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_119
EXIT /B 1
:lbl_119
pushd .
CALL .\check  "%~1"  "schur_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_120
EXIT /B 1
:lbl_120
pushd .
CALL .\check  "%~1"  "mincg_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_121
EXIT /B 1
:lbl_121
pushd .
CALL .\check  "%~1"  "minasa_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_122
EXIT /B 1
:lbl_122
EXIT /B 0
:lbl_63_end
CHDIR _tmp
IF EXIST * DEL /F /Q *
CHDIR ..
IF "%~1"=="msvc" GOTO lbl_123_msvc
GOTO lbl_123___error
:lbl_123_msvc
COPY out\*.h _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_124
echo Error copying ALGLIB library. Did you ran build?
EXIT /B 1
:lbl_124
COPY out\mpn.lib _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_125
echo Error copying ALGLIB library. Did you ran build?
EXIT /B 1
:lbl_125
COPY out\gmp.lib _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_126
echo Error copying ALGLIB library. Did you ran build?
EXIT /B 1
:lbl_126
COPY out\mpfr.lib _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_127
echo Error copying ALGLIB library. Did you ran build?
EXIT /B 1
:lbl_127
COPY out\libmpalglib.lib _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_128
echo Error copying ALGLIB library. Did you ran build?
EXIT /B 1
:lbl_128
GOTO lbl_123___end
:lbl_123___error
echo unknown compiler
EXIT /B 1
:lbl_123___end
IF "%~2"=="ablas" GOTO lbl_129_ablas
IF "%~2"=="ablas_silent" GOTO lbl_129_ablas_silent
IF "%~2"=="ablas_short" GOTO lbl_129_ablas_short
IF "%~2"=="ablasf" GOTO lbl_129_ablasf
IF "%~2"=="ablasf_silent" GOTO lbl_129_ablasf_silent
IF "%~2"=="ablasf_short" GOTO lbl_129_ablasf_short
IF "%~2"=="ortfac" GOTO lbl_129_ortfac
IF "%~2"=="ortfac_silent" GOTO lbl_129_ortfac_silent
IF "%~2"=="ortfac_short" GOTO lbl_129_ortfac_short
IF "%~2"=="hblas" GOTO lbl_129_hblas
IF "%~2"=="hblas_silent" GOTO lbl_129_hblas_silent
IF "%~2"=="hblas_short" GOTO lbl_129_hblas_short
IF "%~2"=="reflections" GOTO lbl_129_reflections
IF "%~2"=="reflections_silent" GOTO lbl_129_reflections_silent
IF "%~2"=="reflections_short" GOTO lbl_129_reflections_short
IF "%~2"=="creflections" GOTO lbl_129_creflections
IF "%~2"=="creflections_silent" GOTO lbl_129_creflections_silent
IF "%~2"=="creflections_short" GOTO lbl_129_creflections_short
IF "%~2"=="sblas" GOTO lbl_129_sblas
IF "%~2"=="sblas_silent" GOTO lbl_129_sblas_silent
IF "%~2"=="sblas_short" GOTO lbl_129_sblas_short
IF "%~2"=="evd" GOTO lbl_129_evd
IF "%~2"=="evd_silent" GOTO lbl_129_evd_silent
IF "%~2"=="evd_short" GOTO lbl_129_evd_short
IF "%~2"=="blas" GOTO lbl_129_blas
IF "%~2"=="blas_silent" GOTO lbl_129_blas_silent
IF "%~2"=="blas_short" GOTO lbl_129_blas_short
IF "%~2"=="rotations" GOTO lbl_129_rotations
IF "%~2"=="rotations_silent" GOTO lbl_129_rotations_silent
IF "%~2"=="rotations_short" GOTO lbl_129_rotations_short
IF "%~2"=="hsschur" GOTO lbl_129_hsschur
IF "%~2"=="hsschur_silent" GOTO lbl_129_hsschur_silent
IF "%~2"=="hsschur_short" GOTO lbl_129_hsschur_short
IF "%~2"=="matgen" GOTO lbl_129_matgen
IF "%~2"=="matgen_silent" GOTO lbl_129_matgen_silent
IF "%~2"=="matgen_short" GOTO lbl_129_matgen_short
IF "%~2"=="hqrnd" GOTO lbl_129_hqrnd
IF "%~2"=="hqrnd_silent" GOTO lbl_129_hqrnd_silent
IF "%~2"=="hqrnd_short" GOTO lbl_129_hqrnd_short
IF "%~2"=="trfac" GOTO lbl_129_trfac
IF "%~2"=="trfac_silent" GOTO lbl_129_trfac_silent
IF "%~2"=="trfac_short" GOTO lbl_129_trfac_short
IF "%~2"=="rcond" GOTO lbl_129_rcond
IF "%~2"=="rcond_silent" GOTO lbl_129_rcond_silent
IF "%~2"=="rcond_short" GOTO lbl_129_rcond_short
IF "%~2"=="trlinsolve" GOTO lbl_129_trlinsolve
IF "%~2"=="trlinsolve_silent" GOTO lbl_129_trlinsolve_silent
IF "%~2"=="trlinsolve_short" GOTO lbl_129_trlinsolve_short
IF "%~2"=="safesolve" GOTO lbl_129_safesolve
IF "%~2"=="safesolve_silent" GOTO lbl_129_safesolve_silent
IF "%~2"=="safesolve_short" GOTO lbl_129_safesolve_short
IF "%~2"=="matinv" GOTO lbl_129_matinv
IF "%~2"=="matinv_silent" GOTO lbl_129_matinv_silent
IF "%~2"=="matinv_short" GOTO lbl_129_matinv_short
IF "%~2"=="bdsvd" GOTO lbl_129_bdsvd
IF "%~2"=="bdsvd_silent" GOTO lbl_129_bdsvd_silent
IF "%~2"=="bdsvd_short" GOTO lbl_129_bdsvd_short
IF "%~2"=="svd" GOTO lbl_129_svd
IF "%~2"=="svd_silent" GOTO lbl_129_svd_silent
IF "%~2"=="svd_short" GOTO lbl_129_svd_short
IF "%~2"=="densesolver" GOTO lbl_129_densesolver
IF "%~2"=="densesolver_silent" GOTO lbl_129_densesolver_silent
IF "%~2"=="densesolver_short" GOTO lbl_129_densesolver_short
IF "%~2"=="xblas" GOTO lbl_129_xblas
IF "%~2"=="xblas_silent" GOTO lbl_129_xblas_silent
IF "%~2"=="xblas_short" GOTO lbl_129_xblas_short
IF "%~2"=="minlbfgs" GOTO lbl_129_minlbfgs
IF "%~2"=="minlbfgs_silent" GOTO lbl_129_minlbfgs_silent
IF "%~2"=="minlbfgs_short" GOTO lbl_129_minlbfgs_short
IF "%~2"=="linmin" GOTO lbl_129_linmin
IF "%~2"=="linmin_silent" GOTO lbl_129_linmin_silent
IF "%~2"=="linmin_short" GOTO lbl_129_linmin_short
IF "%~2"=="odesolver" GOTO lbl_129_odesolver
IF "%~2"=="odesolver_silent" GOTO lbl_129_odesolver_silent
IF "%~2"=="odesolver_short" GOTO lbl_129_odesolver_short
IF "%~2"=="conv" GOTO lbl_129_conv
IF "%~2"=="conv_silent" GOTO lbl_129_conv_silent
IF "%~2"=="conv_short" GOTO lbl_129_conv_short
IF "%~2"=="ftbase" GOTO lbl_129_ftbase
IF "%~2"=="ftbase_silent" GOTO lbl_129_ftbase_silent
IF "%~2"=="ftbase_short" GOTO lbl_129_ftbase_short
IF "%~2"=="fft" GOTO lbl_129_fft
IF "%~2"=="fft_silent" GOTO lbl_129_fft_silent
IF "%~2"=="fft_short" GOTO lbl_129_fft_short
IF "%~2"=="corr" GOTO lbl_129_corr
IF "%~2"=="corr_silent" GOTO lbl_129_corr_silent
IF "%~2"=="corr_short" GOTO lbl_129_corr_short
IF "%~2"=="fht" GOTO lbl_129_fht
IF "%~2"=="fht_silent" GOTO lbl_129_fht_silent
IF "%~2"=="fht_short" GOTO lbl_129_fht_short
IF "%~2"=="autogk" GOTO lbl_129_autogk
IF "%~2"=="autogk_silent" GOTO lbl_129_autogk_silent
IF "%~2"=="autogk_short" GOTO lbl_129_autogk_short
IF "%~2"=="tsort" GOTO lbl_129_tsort
IF "%~2"=="tsort_silent" GOTO lbl_129_tsort_silent
IF "%~2"=="tsort_short" GOTO lbl_129_tsort_short
IF "%~2"=="gammafunc" GOTO lbl_129_gammafunc
IF "%~2"=="gammafunc_silent" GOTO lbl_129_gammafunc_silent
IF "%~2"=="gammafunc_short" GOTO lbl_129_gammafunc_short
IF "%~2"=="gq" GOTO lbl_129_gq
IF "%~2"=="gq_silent" GOTO lbl_129_gq_silent
IF "%~2"=="gq_short" GOTO lbl_129_gq_short
IF "%~2"=="gkq" GOTO lbl_129_gkq
IF "%~2"=="gkq_silent" GOTO lbl_129_gkq_silent
IF "%~2"=="gkq_short" GOTO lbl_129_gkq_short
IF "%~2"=="lsfit" GOTO lbl_129_lsfit
IF "%~2"=="lsfit_silent" GOTO lbl_129_lsfit_silent
IF "%~2"=="lsfit_short" GOTO lbl_129_lsfit_short
IF "%~2"=="minlm" GOTO lbl_129_minlm
IF "%~2"=="minlm_silent" GOTO lbl_129_minlm_silent
IF "%~2"=="minlm_short" GOTO lbl_129_minlm_short
IF "%~2"=="polint" GOTO lbl_129_polint
IF "%~2"=="polint_silent" GOTO lbl_129_polint_silent
IF "%~2"=="polint_short" GOTO lbl_129_polint_short
IF "%~2"=="ratinterpolation" GOTO lbl_129_ratinterpolation
IF "%~2"=="ratinterpolation_silent" GOTO lbl_129_ratinterpolation_silent
IF "%~2"=="ratinterpolation_short" GOTO lbl_129_ratinterpolation_short
IF "%~2"=="ratint" GOTO lbl_129_ratint
IF "%~2"=="ratint_silent" GOTO lbl_129_ratint_silent
IF "%~2"=="ratint_short" GOTO lbl_129_ratint_short
IF "%~2"=="apserv" GOTO lbl_129_apserv
IF "%~2"=="apserv_silent" GOTO lbl_129_apserv_silent
IF "%~2"=="apserv_short" GOTO lbl_129_apserv_short
IF "%~2"=="spline2d" GOTO lbl_129_spline2d
IF "%~2"=="spline2d_silent" GOTO lbl_129_spline2d_silent
IF "%~2"=="spline2d_short" GOTO lbl_129_spline2d_short
IF "%~2"=="spline3" GOTO lbl_129_spline3
IF "%~2"=="spline3_silent" GOTO lbl_129_spline3_silent
IF "%~2"=="spline3_short" GOTO lbl_129_spline3_short
IF "%~2"=="spline1d" GOTO lbl_129_spline1d
IF "%~2"=="spline1d_silent" GOTO lbl_129_spline1d_silent
IF "%~2"=="spline1d_short" GOTO lbl_129_spline1d_short
IF "%~2"=="idwint" GOTO lbl_129_idwint
IF "%~2"=="idwint_silent" GOTO lbl_129_idwint_silent
IF "%~2"=="idwint_short" GOTO lbl_129_idwint_short
IF "%~2"=="nearestneighbor" GOTO lbl_129_nearestneighbor
IF "%~2"=="nearestneighbor_silent" GOTO lbl_129_nearestneighbor_silent
IF "%~2"=="nearestneighbor_short" GOTO lbl_129_nearestneighbor_short
IF "%~2"=="pspline" GOTO lbl_129_pspline
IF "%~2"=="pspline_silent" GOTO lbl_129_pspline_silent
IF "%~2"=="pspline_short" GOTO lbl_129_pspline_short
IF "%~2"=="matdet" GOTO lbl_129_matdet
IF "%~2"=="matdet_silent" GOTO lbl_129_matdet_silent
IF "%~2"=="matdet_short" GOTO lbl_129_matdet_short
IF "%~2"=="sdet" GOTO lbl_129_sdet
IF "%~2"=="sdet_silent" GOTO lbl_129_sdet_silent
IF "%~2"=="sdet_short" GOTO lbl_129_sdet_short
IF "%~2"=="ldlt" GOTO lbl_129_ldlt
IF "%~2"=="ldlt_silent" GOTO lbl_129_ldlt_silent
IF "%~2"=="ldlt_short" GOTO lbl_129_ldlt_short
IF "%~2"=="spdgevd" GOTO lbl_129_spdgevd
IF "%~2"=="spdgevd_silent" GOTO lbl_129_spdgevd_silent
IF "%~2"=="spdgevd_short" GOTO lbl_129_spdgevd_short
IF "%~2"=="sinverse" GOTO lbl_129_sinverse
IF "%~2"=="sinverse_silent" GOTO lbl_129_sinverse_silent
IF "%~2"=="sinverse_short" GOTO lbl_129_sinverse_short
IF "%~2"=="inverseupdate" GOTO lbl_129_inverseupdate
IF "%~2"=="inverseupdate_silent" GOTO lbl_129_inverseupdate_silent
IF "%~2"=="inverseupdate_short" GOTO lbl_129_inverseupdate_short
IF "%~2"=="srcond" GOTO lbl_129_srcond
IF "%~2"=="srcond_silent" GOTO lbl_129_srcond_silent
IF "%~2"=="srcond_short" GOTO lbl_129_srcond_short
IF "%~2"=="ssolve" GOTO lbl_129_ssolve
IF "%~2"=="ssolve_silent" GOTO lbl_129_ssolve_silent
IF "%~2"=="ssolve_short" GOTO lbl_129_ssolve_short
IF "%~2"=="estnorm" GOTO lbl_129_estnorm
IF "%~2"=="estnorm_silent" GOTO lbl_129_estnorm_silent
IF "%~2"=="estnorm_short" GOTO lbl_129_estnorm_short
IF "%~2"=="schur" GOTO lbl_129_schur
IF "%~2"=="schur_silent" GOTO lbl_129_schur_silent
IF "%~2"=="schur_short" GOTO lbl_129_schur_short
IF "%~2"=="mincg" GOTO lbl_129_mincg
IF "%~2"=="mincg_silent" GOTO lbl_129_mincg_silent
IF "%~2"=="mincg_short" GOTO lbl_129_mincg_short
IF "%~2"=="minasa" GOTO lbl_129_minasa
IF "%~2"=="minasa_silent" GOTO lbl_129_minasa_silent
IF "%~2"=="minasa_short" GOTO lbl_129_minasa_short
GOTO lbl_129___error
:lbl_129_ablas
COPY _internal\_run_testablasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testablasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ablas_short
COPY _internal\_run_short_testablasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testablasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ablas_silent
COPY _internal\_run_silent_testablasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testablasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ablasf
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ablasf_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ablasf_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ortfac
COPY _internal\_run_testortfacunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testortfacunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ortfac_short
COPY _internal\_run_short_testortfacunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testortfacunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ortfac_silent
COPY _internal\_run_silent_testortfacunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testortfacunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_hblas
COPY _internal\_run_testhblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testhblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_hblas_short
COPY _internal\_run_short_testhblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testhblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_hblas_silent
COPY _internal\_run_silent_testhblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testhblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_reflections
COPY _internal\_run_testreflectionsunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testreflectionsunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_reflections_short
COPY _internal\_run_short_testreflectionsunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testreflectionsunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_reflections_silent
COPY _internal\_run_silent_testreflectionsunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testreflectionsunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_creflections
COPY _internal\_run_testcreflunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testcreflunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_creflections_short
COPY _internal\_run_short_testcreflunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testcreflunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_creflections_silent
COPY _internal\_run_silent_testcreflunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testcreflunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_sblas
COPY _internal\_run_testsblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_sblas_short
COPY _internal\_run_short_testsblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_sblas_silent
COPY _internal\_run_silent_testsblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_evd
COPY _internal\_run_testevdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testevdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_evd_short
COPY _internal\_run_short_testevdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testevdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_evd_silent
COPY _internal\_run_silent_testevdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testevdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_blas
COPY _internal\_run_testblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_blas_short
COPY _internal\_run_short_testblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_blas_silent
COPY _internal\_run_silent_testblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_rotations
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_rotations_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_rotations_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_hsschur
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_hsschur_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_hsschur_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_matgen
COPY _internal\_run_testmatgenunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmatgenunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_matgen_short
COPY _internal\_run_short_testmatgenunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmatgenunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_matgen_silent
COPY _internal\_run_silent_testmatgenunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmatgenunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_hqrnd
COPY _internal\_run_testhqrndunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testhqrndunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_hqrnd_short
COPY _internal\_run_short_testhqrndunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testhqrndunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_hqrnd_silent
COPY _internal\_run_silent_testhqrndunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testhqrndunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_trfac
COPY _internal\_run_testtrfacunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testtrfacunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_trfac_short
COPY _internal\_run_short_testtrfacunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testtrfacunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_trfac_silent
COPY _internal\_run_silent_testtrfacunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testtrfacunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_rcond
COPY _internal\_run_testrcondunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testrcondunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_rcond_short
COPY _internal\_run_short_testrcondunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testrcondunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_rcond_silent
COPY _internal\_run_silent_testrcondunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testrcondunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_trlinsolve
COPY _internal\_run_testsstunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsstunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_trlinsolve_short
COPY _internal\_run_short_testsstunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsstunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_trlinsolve_silent
COPY _internal\_run_silent_testsstunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsstunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_safesolve
COPY _internal\_run_testsafesolveunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsafesolveunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_safesolve_short
COPY _internal\_run_short_testsafesolveunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsafesolveunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_safesolve_silent
COPY _internal\_run_silent_testsafesolveunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsafesolveunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_matinv
COPY _internal\_run_testmatinvunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmatinvunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_matinv_short
COPY _internal\_run_short_testmatinvunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmatinvunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_matinv_silent
COPY _internal\_run_silent_testmatinvunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmatinvunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_bdsvd
COPY _internal\_run_testbdsvdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testbdsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_bdsvd_short
COPY _internal\_run_short_testbdsvdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testbdsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_bdsvd_silent
COPY _internal\_run_silent_testbdsvdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testbdsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_svd
COPY _internal\_run_testsvdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_svd_short
COPY _internal\_run_short_testsvdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_svd_silent
COPY _internal\_run_silent_testsvdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_densesolver
COPY _internal\_run_testdensesolverunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testdensesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_densesolver_short
COPY _internal\_run_short_testdensesolverunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testdensesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_densesolver_silent
COPY _internal\_run_silent_testdensesolverunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testdensesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_xblas
COPY _internal\_run_testxblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testxblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_xblas_short
COPY _internal\_run_short_testxblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testxblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_xblas_silent
COPY _internal\_run_silent_testxblasunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testxblasunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minlbfgs
COPY _internal\_run_testminlbfgsunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testminlbfgsunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minlbfgs_short
COPY _internal\_run_short_testminlbfgsunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testminlbfgsunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minlbfgs_silent
COPY _internal\_run_silent_testminlbfgsunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testminlbfgsunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_linmin
COPY _internal\_run_testlinminunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testlinminunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_linmin_short
COPY _internal\_run_short_testlinminunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testlinminunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_linmin_silent
COPY _internal\_run_silent_testlinminunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testlinminunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_odesolver
COPY _internal\_run_testodesolverunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testodesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_odesolver_short
COPY _internal\_run_short_testodesolverunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testodesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_odesolver_silent
COPY _internal\_run_silent_testodesolverunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testodesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_conv
COPY _internal\_run_testconvunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testconvunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_conv_short
COPY _internal\_run_short_testconvunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testconvunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_conv_silent
COPY _internal\_run_silent_testconvunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testconvunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ftbase
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ftbase_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ftbase_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_fft
COPY _internal\_run_testfftunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testfftunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_fft_short
COPY _internal\_run_short_testfftunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testfftunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_fft_silent
COPY _internal\_run_silent_testfftunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testfftunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_corr
COPY _internal\_run_testcorrunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testcorrunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_corr_short
COPY _internal\_run_short_testcorrunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testcorrunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_corr_silent
COPY _internal\_run_silent_testcorrunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testcorrunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_fht
COPY _internal\_run_testfhtunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testfhtunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_fht_short
COPY _internal\_run_short_testfhtunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testfhtunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_fht_silent
COPY _internal\_run_silent_testfhtunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testfhtunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_autogk
COPY _internal\_run_testautogk.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testautogk.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_autogk_short
COPY _internal\_run_short_testautogk.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testautogk.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_autogk_silent
COPY _internal\_run_silent_testautogk.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testautogk.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_tsort
COPY _internal\_run_testtsortunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testtsortunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_tsort_short
COPY _internal\_run_short_testtsortunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testtsortunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_tsort_silent
COPY _internal\_run_silent_testtsortunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testtsortunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gammafunc
COPY _internal\_run_testgammaunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgammaunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gammafunc_short
COPY _internal\_run_short_testgammaunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgammaunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gammafunc_silent
COPY _internal\_run_silent_testgammaunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgammaunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gq
COPY _internal\_run_testgq.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgq.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gq_short
COPY _internal\_run_short_testgq.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgq.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gq_silent
COPY _internal\_run_silent_testgq.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgq.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gkq
COPY _internal\_run_testgkq.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgkq.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gkq_short
COPY _internal\_run_short_testgkq.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgkq.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_gkq_silent
COPY _internal\_run_silent_testgkq.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testgkq.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_lsfit
COPY _internal\_run_llstestunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\llstestunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_lsfit_short
COPY _internal\_run_short_llstestunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\llstestunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_lsfit_silent
COPY _internal\_run_silent_llstestunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\llstestunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minlm
COPY _internal\_run_testminlmunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testminlmunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minlm_short
COPY _internal\_run_short_testminlmunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testminlmunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minlm_silent
COPY _internal\_run_silent_testminlmunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testminlmunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_polint
COPY _internal\_run_testpolintunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testpolintunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_polint_short
COPY _internal\_run_short_testpolintunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testpolintunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_polint_silent
COPY _internal\_run_silent_testpolintunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testpolintunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ratinterpolation
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ratinterpolation_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ratinterpolation_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ratint
COPY _internal\_run_testratinterpolation.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testratinterpolation.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ratint_short
COPY _internal\_run_short_testratinterpolation.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testratinterpolation.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ratint_silent
COPY _internal\_run_silent_testratinterpolation.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testratinterpolation.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_apserv
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_apserv_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_apserv_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_spline2d
COPY _internal\_run_testspline2dunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspline2dunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spline2d_short
COPY _internal\_run_short_testspline2dunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspline2dunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spline2d_silent
COPY _internal\_run_silent_testspline2dunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspline2dunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spline3
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_spline3_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_spline3_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_spline1d
COPY _internal\_run_testspline1dunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspline1dunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spline1d_short
COPY _internal\_run_short_testspline1dunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspline1dunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spline1d_silent
COPY _internal\_run_silent_testspline1dunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspline1dunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_idwint
COPY _internal\_run_testidwunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testidwunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_idwint_short
COPY _internal\_run_short_testidwunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testidwunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_idwint_silent
COPY _internal\_run_silent_testidwunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testidwunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_nearestneighbor
COPY _internal\_run_testnearestneighborunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testnearestneighborunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_nearestneighbor_short
COPY _internal\_run_short_testnearestneighborunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testnearestneighborunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_nearestneighbor_silent
COPY _internal\_run_silent_testnearestneighborunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testnearestneighborunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_pspline
COPY _internal\_run_testpsplineunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testpsplineunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_pspline_short
COPY _internal\_run_short_testpsplineunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testpsplineunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_pspline_silent
COPY _internal\_run_silent_testpsplineunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testpsplineunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_matdet
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_matdet_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_matdet_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_sdet
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_sdet_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_sdet_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ldlt
COPY _internal\_run_testldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ldlt_short
COPY _internal\_run_short_testldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ldlt_silent
COPY _internal\_run_silent_testldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spdgevd
COPY _internal\_run_testspdgevdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspdgevdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spdgevd_short
COPY _internal\_run_short_testspdgevdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspdgevdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_spdgevd_silent
COPY _internal\_run_silent_testspdgevdunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testspdgevdunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_sinverse
COPY _internal\_run_testinvldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testinvldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_sinverse_short
COPY _internal\_run_short_testinvldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testinvldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_sinverse_silent
COPY _internal\_run_silent_testinvldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testinvldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_inverseupdate
COPY _internal\_run_testinverseupdateunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testinverseupdateunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_inverseupdate_short
COPY _internal\_run_short_testinverseupdateunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testinverseupdateunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_inverseupdate_silent
COPY _internal\_run_silent_testinverseupdateunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testinverseupdateunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_srcond
COPY _internal\_run_testrcondldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testrcondldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_srcond_short
COPY _internal\_run_short_testrcondldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testrcondldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_srcond_silent
COPY _internal\_run_silent_testrcondldltunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testrcondldltunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_ssolve
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ssolve_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_ssolve_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_estnorm
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_129___end
:lbl_129_estnorm_short
EXIT /B 0
GOTO lbl_129___end
:lbl_129_estnorm_silent
EXIT /B 0
GOTO lbl_129___end
:lbl_129_schur
COPY _internal\_run_testschurunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testschurunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_schur_short
COPY _internal\_run_short_testschurunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testschurunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_schur_silent
COPY _internal\_run_silent_testschurunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testschurunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_mincg
COPY _internal\_run_testmincgunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmincgunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_mincg_short
COPY _internal\_run_short_testmincgunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmincgunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_mincg_silent
COPY _internal\_run_silent_testmincgunit.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testmincgunit.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minasa
COPY _internal\_run_testasa.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testasa.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minasa_short
COPY _internal\_run_short_testasa.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testasa.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129_minasa_silent
COPY _internal\_run_silent_testasa.cpp _tmp\_test.cpp > NUL 2> NUL
COPY tests\testasa.* _tmp > NUL 2> NUL
GOTO lbl_129___end
:lbl_129___error
echo unknown unit
EXIT /B 1
:lbl_129___end
IF "%~1"=="msvc" GOTO lbl_130_msvc
GOTO lbl_130___error
:lbl_130_msvc
CHDIR _tmp
cl /nologo /EHsc /I. /Fe_test.exe %~3 *.cpp  libmpalglib.lib mpfr.lib gmp.lib >> ../log.txt 2>&1
IF NOT ERRORLEVEL 1 GOTO lbl_131
echo Error while compiling (see ../log.txt for more info)
CHDIR ..
EXIT /B 1
:lbl_131
CHDIR ..
pushd _tmp
.\_test
popd
GOTO lbl_130___end
:lbl_130___error
echo unknown compiler
EXIT /B 1
:lbl_130___end
EXIT /B 0
