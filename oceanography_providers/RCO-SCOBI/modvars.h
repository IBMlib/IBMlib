
!
!       RCO-SCOBI
!
!---------------------------------------------------------------------
!
! variable descriptions, units and grid indicator (see varinfo.h)
! This is used by genmoddata to generate the moddata.h file
! Do not remove or add lines in the middle of the list!
! Maybe this data should be moved to varinfo.h ?

!     This file should be formated so it works well to include in both
!     fixed format (aka F77) and free format (aka F90) Fortran files.
!     At the moment it means some care must be used for continuation
!     lines and comments.
!
!     Comments can start anywhere on a line but must be preceeded by an
!     exclamaition mark (!) even if they start at or after position 72.
!
!     Lines to be continued must have a ampersand (&) in position 72 or
!     later (so it is regarded as comment for fixed format and not a
!     comment in free format) and the continuation line must have an
!     ampersand (&) in position 6.

      data (varname(n),varunit(n),vartgrd(n),n=1,99)/                   &
     &  'free surface height      ', 'cm          ', .true. ,           &! 1
     &  'barotropic u-velocity    ', 'cm/s        ', .false. ,          &! 2
     &  'barotropic v-velocity    ', 'cm/s        ', .false. ,          &! 3
     &  'potential temperature    ', 'C           ', .true. ,           &! 4
     &  'salinity                 ', 'psu         ', .true. ,           &! 5
     &  'turbulent kinetic energy ', 'cm**2/s**2  ', .true. ,           &! 6
     &  'dissipation              ', 'cm**2/s**3  ', .true. ,           &! 7
     &  'u-velocity               ', 'cm/s        ', .false. ,          &! 8
     &  'v-velocity               ', 'cm/s        ', .false. ,          &! 9
     &  'zooplankton              ', 'mgC/m**3    ', .true. ,           &! 10
     &  'phytoplankton 1          ', 'mgCHL/m**3  ', .true. ,           &! 11
     &  'phytoplankton 2          ', 'mgCHL/m**3  ', .true. ,           &! 12
     &  'phytoplankton 3          ', 'mgCHL/m**3  ', .true. ,           &! 13
     &  'detritus                 ', 'mgC/m**3    ', .true. ,           &! 14
     &  'ammonimum                ', 'mmolN/m**3  ', .true. ,           &! 15
     &  'nitrate                  ', 'mmolN/m**3  ', .true. ,           &! 16
     &  'phosphate                ', 'mmolP/m**3  ', .true. ,           &! 17
     &  'oxygen                   ', 'ml/l        ', .true. ,           &! 18
     &  'benthic nitrogen         ', 'mmolN/m**2  ', .true. ,           &! 19
     &  'benthic phosphorus       ', 'mmolP/m**2  ', .true. ,           &! 20
     &  'pelagic resusp tracer    ', 'mg C/m**3   ', .true. ,           &! 21
     &  'benthic resusp tracer    ', 'mmolN/m**2  ', .true. ,           &! 22
     &  'pelagic IronPO4          ', 'mmolP/m**3  ', .true. ,           &! 23
     &  'benthic IronPO4          ', 'mmolP/m**2  ', .true. ,           &! 24
     &  'phosphorus detritus      ', 'mgC/m**3    ', .true. ,           &! 25
     &  'pelagic akinetes         ', 'mmolN/m**3  ', .true. ,           &! 26
     &  'pelagic recuting cells   ', 'mmolN/m**3  ', .true. ,           &! 27
     &  'benthic akinetes         ', 'mmolN/m**2  ', .true. ,           &! 28
     &  'benthic maturation time  ', 'diagnostic  ', .true. ,           &! 29
     &  'clc mean growth rate     ', '1/d         ', .true. ,           &! 30
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &! 31
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &! 32
     &  'htop acc                 ', '            ', .true. ,           &! 33
     &  'hbot acc                 ', '            ', .true. ,           &! 34
     &  'old free surface height  ', 'cm          ', .true. ,           &! 35
     &  'Sigma11a ucomp triangle1 ', '  dyn/cm    ', .true. ,           &! 36
     &  'Sigma11b ucomp triangle2 ', '  dyn/cm    ', .true. ,           &! 37
     &  'Sigma11c ucomp triangle3 ', '  dyn/cm    ', .true. ,           &! 38
     &  'Sigma11d ucomp triangle4 ', '  dyn/cm    ', .true. ,           &! 39
     &  'Sigma12a shear triangle1 ', '  dyn/cm    ', .true. ,           &! 40
     &  'Sigma12b shear triangle2 ', '  dyn/cm    ', .true. ,           &! 41
     &  'Sigma12c shear triangle3 ', '  dyn/cm    ', .true. ,           &! 42
     &  'Sigma12d shear triangle4 ', '  dyn/cm    ', .true. ,           &! 43
     &  'Sigma22a vcomp triangle1 ', '  dyn/cm    ', .true. ,           &! 44
     &  'Sigma22b vcomp triangle2 ', '  dyn/cm    ', .true. ,           &! 45
     &  'Sigma22c vcomp triangle3 ', '  dyn/cm    ', .true. ,           &! 46
     &  'Sigma22d vcomp triangle4 ', '  dyn/cm    ', .true. ,           &! 47
     &  'heat content of brine    ', '  null      ', .true. ,           &! 48
     &  'number of ice layers     ', '            ', .true. ,           &! 49
     &  'fraction of ice area     ', '  null      ', .true. ,           &! 50
     &  'thickness of ice layer   ', 'cm          ', .true. ,           &! 51
     &  'surface temp of ice      ', 'C           ', .true. ,           &! 52
     &  'temp of first ice layer  ', 'C           ', .true. ,           &! 53
     &  'smf x                    ', 'dyn/cm**2   ', .false. ,          &! 54
     &  'smf y                    ', 'dyn/cm**2   ', .false. ,          &! 55
     &  'uice x                   ', 'cm/s        ', .false. ,          &! 56
     &  'vice y                   ', 'cm/s        ', .false. ,          &! 57
     &  'hsno                     ', '            ', .true. ,           &! 58
     &  'tsno                     ', '            ', .true. ,           &! 59
     &  'tice2                    ', '            ', .true. ,           &! 60
     &  'volice                   ', '            ', .true. ,           &! 61
     &  'htop                     ', '            ', .true. ,           &! 62
     &  'hbot                     ', '            ', .true. ,           &! 63
     &  'qs_a_w_A                 ', '            ', .true. ,           &! 64
     &  'qbin_a_w_A               ', '            ', .true. ,           &! 65
     &  'qbout_a_w_A              ', '            ', .true. ,           &! 66
     &  'qk_a_w_A                 ', '            ', .true. ,           &! 67
     &  'qe_a_w_A                 ', '            ', .true. ,           &! 68
     &  'qs_a_i_A                 ', '            ', .true. ,           &! 69
     &  'qbin_a_i_A               ', '            ', .true. ,           &! 70
     &  'qbout_a_i_A              ', '            ', .true. ,           &! 71
     &  'qk_a_i_A                 ', '            ', .true. ,           &! 72
     &  'qe_a_i_A                 ', '            ', .true. ,           &! 73
     &  'qs_w_i_O                 ', '            ', .true. ,           &! 74
     &  'qbot_w_i_O               ', '            ', .true. ,           &! 75
     &  'qs_a_w_O                 ', '            ', .true. ,           &! 76
     &  'qbin_a_w_O               ', '            ', .true. ,           &! 77
     &  'qbout_a_w_O              ', '            ', .true. ,           &! 78
     &  'qk_a_w_O                 ', '            ', .true. ,           &! 79
     &  'qe_a_w_O                 ', '            ', .true. ,           &! 80
     &  'qs_a_i                   ', '            ', .true. ,           &! 81
     &  'qbin_a_i                 ', '            ', .true. ,           &! 82
     &  'qbout_a_i                ', '            ', .true. ,           &! 83
     &  'qk_a_i                   ', '            ', .true. ,           &! 84
     &  'qe_a_i                   ', '            ', .true. ,           &! 85
     &  'qs_w_i                   ', '            ', .true. ,           &! 86
     &  'qbot_w_i                 ', '            ', .true. ,           &! 87
     &  'volice_freez             ', '            ', .true. ,           &! 88
     &  'volice_melt              ', '            ', .true. ,           &! 89
     &  'htop_melt                ', '            ', .true. ,           &! 90
     &  'hbot_freez               ', '            ', .true. ,           &! 91
     &  'hbot_melt                ', '            ', .true. ,           &! 92
     &  'pme2d                    ', '            ', .true. ,           &! 93
     &  'precip                   ', '            ', .true. ,           &! 94
     &  'runoff                   ', '            ', .true. ,           &! 95
     &  'qs_tot_wi                ', '            ', .true. ,           &! 96
     &  'salinity_gr              ', '            ', .true. ,           &! 97
     &  'volume_gr                ', '            ', .true. ,           &! 98
     &  'rafted ice thickness     ', 'cm          ', .true. /            ! 99
      data (varname(n),varunit(n),vartgrd(n),n=100,198)/                &
     &  'rafted ice concentration ', 'fraction    ', .true. ,           &!100
     &  'ridged ice thickness     ', 'cm          ', .true. ,           &!101
     &  'ridged ice concentration ', 'fraction    ', .true. ,           &!102
     &  'slp                      ', '  Pa -100000', .true. ,           &!103
     &  'u10                      ', '  m/s       ', .true. ,           &!104
     &  'v10                      ', '  m/s       ', .true. ,           &!105
     &  'tair                     ', '  Deg.C     ', .true. ,           &!106
     &  'rqa                      ', '            ', .true. ,           &!107
     &  'cloudiness               ', '            ', .true. ,           &!108
     &  'sumvar1                  ', '            ', .true. ,           &!109
     &  'sumvar2                  ', '            ', .true. ,           &!110
     &  'sumvar3                  ', '            ', .true. ,           &!111
     &  'sumvar4                  ', '            ', .true. ,           &!112
     &  'sumvar5                  ', '            ', .true. ,           &!113
     &  'sumvar6                  ', '            ', .true. ,           &!114
     &  'sumvar7,nr resusp events ', '            ', .true. ,           &!115
     &  'sumvar8,max TAU          ', '  N/m2      ', .true. ,           &!116
     &  'sumvar9,95 % tau         ', '            ', .true. ,           &!117
     &  'sumvar10,sum of all tau  ', '  N/m2      ', .true. ,           &!118
     &  'sumvar11,average tau     ', '  N/m2      ', .true. ,           &!119
     &  'sumvar12,99 % tau        ', '            ', .true. ,           &!120
     &  'sumvar13,tau due to waves', '  N/m2      ', .true. ,           &!121
     &  'sumvar14,tau due to curr ', '  N/m2      ', .true. ,           &!122
     &  'sumvar15,tau max         ', '  N/m2      ', .true. ,           &!123
     &  'tausum1                  ', '  nr        ', .true. ,           &!124
     &  'tausum2                  ', '  nr        ', .true. ,           &!125
     &  'tausum3                  ', '  nr        ', .true. ,           &!126
     &  'tausum4                  ', '  nr        ', .true. ,           &!127
     &  'tausum5                  ', '  nr        ', .true. ,           &!128
     &  'tausum6                  ', '  nr        ', .true. ,           &!129
     &  'tausum7                  ', '  nr        ', .true. ,           &!130
     &  'tausum8                  ', '  nr        ', .true. ,           &!131
     &  'tausum9                  ', '  nr        ', .true. ,           &!132
     &  'tausum10                 ', '  nr        ', .true. ,           &!133
     &  'tausum11                 ', '  nr        ', .true. ,           &!134
     &  'tausum12                 ', '  nr        ', .true. ,           &!135
     &  'tausum13                 ', '  nr        ', .true. ,           &!136
     &  'tausum14                 ', '  nr        ', .true. ,           &!137
     &  'tausum15                 ', '  nr        ', .true. ,           &!138
     &  'tausum16                 ', '  nr        ', .true. ,           &!139
     &  'tausum17                 ', '  nr        ', .true. ,           &!140
     &  'tausum18                 ', '  nr        ', .true. ,           &!141
     &  'tausum19                 ', '  nr        ', .true. ,           &!142
     &  'tausum20                 ', '  nr        ', .true. ,           &!143
     &  'tausum21                 ', '  nr        ', .true. ,           &!144
     &  'tausum22                 ', '  nr        ', .true. ,           &!145
     &  'tausum23                 ', '  nr        ', .true. ,           &!146
     &  'tausum24                 ', '  nr        ', .true. ,           &!147
     &  'tausum25                 ', '  nr        ', .true. ,           &!148
     &  'tausum26                 ', '  nr        ', .true. ,           &!149
     &  'tausum27                 ', '  nr        ', .true. ,           &!150
     &  'tausum28                 ', '  nr        ', .true. ,           &!151
     &  'tausum29                 ', '  nr        ', .true. ,           &!152
     &  'tausum30                 ', '  nr        ', .true. ,           &!153
     &  'tausum31                 ', '  nr        ', .true. ,           &!154
     &  'tausum32                 ', '  nr        ', .true. ,           &!155
     &  'tausum33                 ', '  nr        ', .true. ,           &!156
     &  'tausum34                 ', '  nr        ', .true. ,           &!157
     &  'tausum35                 ', '  nr        ', .true. ,           &!158
     &  'tausum36                 ', '  nr        ', .true. ,           &!159
     &  'tausum37                 ', '  nr        ', .true. ,           &!160
     &  'tausum38                 ', '  nr        ', .true. ,           &!161
     &  'tausum39                 ', '  nr        ', .true. ,           &!162
     &  'tausum40                 ', '  nr        ', .true. ,           &!163
     &  'tausum41                 ', '  nr        ', .true. ,           &!164
     &  'hwave signif. wave height', '  cm        ', .true. ,           &!165
     &  'sumvar16, secchi scobi   ', '  m         ', .true. ,           &!166
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!167
     &  'Primary production       ', '  kton C    ', .true. ,           &!168
     &  'Nitrogen fixation        ', '  kton N    ', .true. ,           &!169
     &  'Pel DIN uptake to pp     ', '  kton N    ', .true. ,           &!170
     &  'Pel DIP uptake to pp     ', '  kton P    ', .true. ,           &!171
     &  'Pel N remineralization   ', '  kton N    ', .true. ,           &!172
     &  'Pel P remineralization   ', '  kton P    ', .true. ,           &!173
     &  'Net sedimentation of N   ', '  kton N    ', .true. ,           &!174
     &  'Net sedimentation of P   ', '  kton P    ', .true. ,           &!175
     &  'Pelagic denitrification  ', '  kton N    ', .true. ,           &!176
     &  'Benthic denitrification  ', '  kton P    ', .true. ,           &!177
     &  'Benthic N release        ', '  kton N    ', .true. ,           &!178
     &  'Benthic P release        ', '  kton P    ', .true. ,           &!179
     &  'Burial of N              ', '  kton N    ', .true. ,           &!180
     &  'Burial of P              ', '  kton P    ', .true. ,           &!181
     &  'Atmospheric supply P     ', '  kton P    ', .true. ,           &!182
     &  'Atmospheric supply N     ', '  kton N    ', .true. ,           &!183
     &  'River load bioavailable P', '  kton P    ', .true. ,           &!184
     &  'River load bioavailable N', '  kton N    ', .true. ,           &!185
     &  'Point load bioavailable P', '  kton P    ', .true. ,           &!186
     &  'Point load bioavailable N', '  kton N    ', .true. ,           &!187
     &  'River runoff             ', '  km3       ', .true. ,           &!188
     &  'River load DIP           ', '  kton P    ', .true. ,           &!189
     &  'River load DIN           ', '  kton N    ', .true. ,           &!190
     &  'Point load DIP           ', '  kton P    ', .true. ,           &!191
     &  'Point load DIN           ', '  kton N    ', .true. ,           &!192
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!193
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!194
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!195
     &  'ice stress div comp s11  ', 'dyn/cm^2    ', .true. ,           &!196
     &  'ice stress div comp s12  ', 'dyn/cm^2    ', .true. ,           &!197
     &  'ice stress div comp s21  ', 'dyn/cm^2    ', .true. /            !198
      data (varname(n),varunit(n),vartgrd(n),n=199,400)/                &
     &  'ice stress div comp s22  ', 'dyn/cm^2    ', .true. ,           &!199
     &  'level ice1 conc          ', 'fraction    ', .true. ,           &!200
     &  'level ice2 conc          ', 'fraction    ', .true. ,           &!201
     &  'level ice3 conc          ', 'fraction    ', .true. ,           &!202
     &  'level ice4 conc          ', 'fraction    ', .true. ,           &!203
     &  'level ice5 conc          ', 'fraction    ', .true. ,           &!204
     &  'level ice6 conc          ', 'fraction    ', .true. ,           &!205
     &  'level ice7 conc          ', 'fraction    ', .true. ,           &!206
     &  'level ice8 conc          ', 'fraction    ', .true. ,           &!207
     &  'level ice9 conc          ', 'fraction    ', .true. ,           &!208
     &  'level ice10 conc         ', 'fraction    ', .true. ,           &!209
     &  'ridged ice mass          ', '?           ', .true. ,           &!210
     &  'rafted ice mass          ', '?           ', .true. ,           &!211
     &  'level ice1 mean thickness', 'cm          ', .true. ,           &!212
     &  'level ice2 mean thickness', 'cm          ', .true. ,           &!213
     &  'level ice3 mean thickness', 'cm          ', .true. ,           &!214
     &  'level ice4 mean thickness', 'cm          ', .true. ,           &!215
     &  'level ice5 mean thickness', 'cm          ', .true. ,           &!216
     &  'level ice6 mean thickness', 'cm          ', .true. ,           &!217
     &  'level ice7 mean thickness', 'cm          ', .true. ,           &!218
     &  'level ice8 mean thickness', 'cm          ', .true. ,           &!219
     &  'level ice9 mean thickness', 'cm          ', .true. ,           &!220
     &  'level ice10 meanThickness', 'cm          ', .true. ,           &!221
     &  'ridged snow mass         ', '            ', .true. ,           &!222
     &  'rafted snow mass         ', '            ', .true. ,           &!223
     &  'level snow1 mass         ', '            ', .true. ,           &!224
     &  'level snow2 mass         ', '            ', .true. ,           &!225
     &  'level snow3 mass         ', '            ', .true. ,           &!226
     &  'level snow4 mass         ', '            ', .true. ,           &!227
     &  'level snow5 mass         ', '            ', .true. ,           &!228
     &  'level snow6 mass         ', '            ', .true. ,           &!229
     &  'level snow7 mass         ', '            ', .true. ,           &!230
     &  'level snow8 mass         ', '            ', .true. ,           &!231
     &  'level snow9 mass         ', '            ', .true. ,           &!232
     &  'level snow10 mass        ', '            ', .true. ,           &!233
     &  'ridged heat brine        ', '  null      ', .true. ,           &!234
     &  'rafted heat brine        ', '  null      ', .true. ,           &!235
     &  'level ice1 heat brine    ', '  null      ', .true. ,           &!236
     &  'level ice2 heat brine    ', '  null      ', .true. ,           &!237
     &  'level ice3 heat brine    ', '  null      ', .true. ,           &!238
     &  'level ice4 heat brine    ', '  null      ', .true. ,           &!239
     &  'level ice5 heat brine    ', '  null      ', .true. ,           &!240
     &  'level ice6 heat brine    ', '  null      ', .true. ,           &!241
     &  'level ice7 heat brine    ', '  null      ', .true. ,           &!242
     &  'level ice8 heat brine    ', '  null      ', .true. ,           &!243
     &  'level ice9 heat brine    ', '  null      ', .true. ,           &!244
     &  'level ice10 heat brine   ', '  null      ', .true. ,           &!245
     &  'ridged surface temp      ', '  null      ', .true. ,           &!246
     &  'rafted surface temp      ', '  null      ', .true. ,           &!247
     &  'level ice1 surface temp  ', '  null      ', .true. ,           &!248
     &  'level ice2 surface temp  ', '  null      ', .true. ,           &!249
     &  'level ice3 surface temp  ', '  null      ', .true. ,           &!250
     &  'level ice4 surface temp  ', '  null      ', .true. ,           &!251
     &  'level ice5 surface temp  ', '  null      ', .true. ,           &!252
     &  'level ice6 surface temp  ', '  null      ', .true. ,           &!253
     &  'level ice7 surface temp  ', '  null      ', .true. ,           &!254
     &  'level ice8 surface temp  ', '  null      ', .true. ,           &!255
     &  'level ice9 surface temp  ', '  null      ', .true. ,           &!256
     &  'level ice10 surface temp ', '  null      ', .true. ,           &!257
     &  'ridged snow temp         ', '  null      ', .true. ,           &!258
     &  'rafted snow temp         ', '  null      ', .true. ,           &!259
     &  'level ice 1 snow temp    ', '  null      ', .true. ,           &!260
     &  'level ice 2 snow temp    ', '  null      ', .true. ,           &!261
     &  'level ice 3 snow temp    ', '  null      ', .true. ,           &!262
     &  'level ice 4 snow temp    ', '  null      ', .true. ,           &!263
     &  'level ice 5 snow temp    ', '  null      ', .true. ,           &!264
     &  'level ice 6 snow temp    ', '  null      ', .true. ,           &!265
     &  'level ice 7 snow temp    ', '  null      ', .true. ,           &!266
     &  'level ice 8 snow temp    ', '  null      ', .true. ,           &!267
     &  'level ice 9 snow temp    ', '  null      ', .true. ,           &!268
     &  'level ice 10 snow temp   ', '  null      ', .true. ,           &!269
     &  'ridged ice temp 1        ', '  null      ', .true. ,           &!270
     &  'ridged ice temp 2        ', '  null      ', .true. ,           &!271
     &  'rafted ice temp 1        ', '  null      ', .true. ,           &!272
     &  'rafted ice temp 2        ', '  null      ', .true. ,           &!273
     &  'level ice1 temp 1        ', '  null      ', .true. ,           &!274
     &  'level ice1 temp 2        ', '  null      ', .true. ,           &!275
     &  'level ice2 temp 1        ', '  null      ', .true. ,           &!276
     &  'level ice2 temp 2        ', '  null      ', .true. ,           &!277
     &  'level ice3 temp 1        ', '  null      ', .true. ,           &!278
     &  'level ice3 temp 2        ', '  null      ', .true. ,           &!279
     &  'level ice4 temp 1        ', '  null      ', .true. ,           &!280
     &  'level ice4 temp 2        ', '  null      ', .true. ,           &!281
     &  'level ice5 temp 1        ', '  null      ', .true. ,           &!282
     &  'level ice5 temp 2        ', '  null      ', .true. ,           &!283
     &  'level ice6 temp 1        ', '  null      ', .true. ,           &!284
     &  'level ice6 temp 2        ', '  null      ', .true. ,           &!285
     &  'level ice7 temp 1        ', '  null      ', .true. ,           &!286
     &  'level ice7 temp 2        ', '  null      ', .true. ,           &!287
     &  'level ice8 temp 1        ', '  null      ', .true. ,           &!288
     &  'level ice8 temp 2        ', '  null      ', .true. ,           &!289
     &  'level ice9 temp 1        ', '  null      ', .true. ,           &!290
     &  'level ice9 temp 2        ', '  null      ', .true. ,           &!291
     &  'level ice10 temp 1       ', '  null      ', .true. ,           &!292
     &  'level ice10 temp 2       ', '  null      ', .true. ,           &!293
     &  'ridged number of layers  ', '            ', .true. ,           &!294
     &  'rafted number of layers  ', '            ', .true. ,           &!295
     &  'ice1 number of layers    ', '            ', .true. ,           &!296
     &  'ice2 number of layers    ', '            ', .true. ,           &!297
     &  'ice3 number of layers    ', '            ', .true. ,           &!298
     &  'ice4 number of layers    ', '            ', .true. ,           &!299
     &  'ice5 number of layers    ', '            ', .true. ,           &!300
     &  'ice6 number of layers    ', '            ', .true. ,           &!301
     &  'ice7 number of layers    ', '            ', .true. ,           &!302
     &  'ice8 number of layers    ', '            ', .true. ,           &!303
     &  'ice9 number of layers    ', '            ', .true. ,           &!304
     &  'ice10 number of layers   ', '            ', .true. ,           &!305
     &  'uice sum                 ', 'm/s         ', .false.,           &!306
     &  'vice sum                 ', 'm/s         ', .false.,           &!307
     &  'stf1_2d                  ', 'K/s         ', .true. ,           &!308
     &  'stf2_2d                  ', 'S/s         ', .true. ,           &!309
     &  'ice strength             ', 'bla         ', .true. ,           &!310
     &  'albedo snow + ice        ', '%           ', .true. ,           &!311
     &  'albedois                 ', '%           ', .true. ,           &!312
     &  'transport org P north    ', 'ngP         ', .true. ,           &!313
     &  'transport org P south    ', 'ngP         ', .true. ,           &!314
     &  'transport org P east     ', 'ngP         ', .true. ,           &!315
     &  'transport org P west     ', 'ngP         ', .true. ,           &!316
     &  'transport din north      ', 'nmolN       ', .true. ,           &!317
     &  'transport din south      ', 'nmolN       ', .true. ,           &!318
     &  'transport din east       ', 'nmolN       ', .true. ,           &!319
     &  'transport din west       ', 'nmolN       ', .true. ,           &!320
     &  'transport dip north      ', 'nmolP       ', .true. ,           &!321
     &  'transport dip south      ', 'nmolP       ', .true. ,           &!322
     &  'transport dip east       ', 'nmolP       ', .true. ,           &!323
     &  'transport dip west       ', 'nmolP       ', .true. ,           &!324
     &  'transp pelag resusp north', 'ngC         ', .true. ,           &!325
     &  'transp pelag resusp south', 'ngC         ', .true. ,           &!326
     &  'transp pelag resusp east ', 'ngC         ', .true. ,           &!327
     &  'transp pelag resusp west ', 'ngC         ', .true. ,           &!328
     &  'divergence ice stress x  ', '  dyn/cm^2  ', .true. ,           &!329
     &  'divergence ice stress y  ', '  dyn/cm^2  ', .true. ,           &!330
     &  'divergence ice vel.      ', '  null      ', .true. ,           &!331
     &  'transport ironp north    ', 'nmolP       ', .true. ,           &!332
     &  'transport ironp south    ', 'nmolP       ', .true. ,           &!333
     &  'transport ironp east     ', 'nmolP       ', .true. ,           &!334
     &  'transport ironp west     ', 'nmolP       ', .true. ,           &!335
     &  'transport pdt north      ', 'nmolP       ', .true. ,           &!336
     &  'transport pdt south      ', 'nmolP       ', .true. ,           &!337
     &  'transport pdt east       ', 'nmolP       ', .true. ,           &!338
     &  'transport pdt west       ', 'nmolP       ', .true. ,           &!339
     &  'clc mattime              ', 'days        ', .true. ,           &!340
     &  'clc akimax               ', 'mmolN/m**2  ', .true. ,           &!341
     &  'clc akisum               ', 'mmolN/m**2  ', .true. ,           &!342
     &  'OXYGEN PRODUCTION        ', 'kton oxy    ', .true. ,           &!343
     &  'PEL OXY CONSUMPT         ', 'kton oxy    ', .true. ,           &!344
     &  'BEN OXY CONSUMPT         ', 'kton oxy    ', .true. ,           &!345
     &  'PEL SO4 RED (neg oxy)    ', 'kton oxy    ', .true. ,           &!346
     &  'BEN SO4 RED (neg oxy)    ', 'kton oxy    ', .true. ,           &!347
     &  'DIATOM PRODUCTION        ', 'kton C      ', .true. ,           &!348
     &  'FLA & al PRODUCTION      ', 'kton C      ', .true. ,           &!349
     &  'CYANO PRODUCTION         ', 'kton C      ', .true. ,           &!350
     &  'DIATOM GRAZINGE          ', 'kton C      ', .true. ,           &!351
     &  'FLA & al GRAZING         ', 'kton C      ', .true. ,           &!352
     &  'CYA GRAZING              ', 'kton C      ', .true. ,           &!353
     &  'DET GRAZING              ', 'kton C      ', .true. ,           &!354
     &  'PEDATION ON ZOOPLANKTON  ', 'kton C      ', .true. ,           &!355
     &  'ZOO NH4 EXCRETION        ', 'kton C      ', .true. ,           &!356
     &  'ZOOPLANKTON FEACES       ', 'kton C      ', .true. ,           &!357
     &  'w-velocity               ', 'cm/s        ', .true. ,           &!358
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!359
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!360
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!361
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!362
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!363
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!364
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!365
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!366
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!367
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!368
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!369
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!370
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!371
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!372
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!373
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!374
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!375
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!376
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!377
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!378
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!379
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!380
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!381
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!382
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!383
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!384
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!385
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!386
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!387
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!388
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!389
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!390
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!391
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!392
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!393
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!394
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!395
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!396
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!397
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!398
     &  'UNUSED VARIABLE          ', '  null      ', .true. ,           &!399
     &  'UNUSED VARIABLE          ', '  null      ', .true. /            !400

      data (varname(n),varunit(n),vartgrd(n),n=401,500)/                &
     &  'tracer 1                 ', '            ', .true.,            &!401
     &  'tracer 2                 ', '            ', .true.,            &!402
     &  'tracer 3                 ', '            ', .true.,            &!403
     &  'tracer 4                 ', '            ', .true.,            &!404
     &  'tracer 5                 ', '            ', .true.,            &!405
     &  'tracer 6                 ', '            ', .true.,            &!406
     &  'tracer 7                 ', '            ', .true.,            &!407
     &  'tracer 8                 ', '            ', .true.,            &!408
     &  'tracer 9                 ', '            ', .true.,            &!409
     &  'tracer 10                ', '            ', .true.,            &!410
     &  'tracer 11                ', '            ', .true.,            &!411
     &  'tracer 12                ', '            ', .true.,            &!412
     &  'tracer 13                ', '            ', .true.,            &!413
     &  'tracer 14                ', '            ', .true.,            &!414
     &  'tracer 15                ', '            ', .true.,            &!415
     &  'tracer 16                ', '            ', .true.,            &!416
     &  'tracer 17                ', '            ', .true.,            &!417
     &  'tracer 18                ', '            ', .true.,            &!418
     &  'tracer 19                ', '            ', .true.,            &!419
     &  'tracer 20                ', '            ', .true.,            &!420
     &  'tracer 21                ', '            ', .true.,            &!421
     &  'tracer 22                ', '            ', .true.,            &!422
     &  'tracer 23                ', '            ', .true.,            &!423
     &  'tracer 24                ', '            ', .true.,            &!424
     &  'tracer 25                ', '            ', .true.,            &!425
     &  'tracer 26                ', '            ', .true.,            &!426
     &  'tracer 27                ', '            ', .true.,            &!427
     &  'tracer 28                ', '            ', .true.,            &!428
     &  'tracer 29                ', '            ', .true.,            &!429
     &  'tracer 30                ', '            ', .true.,            &!430
     &  'tracer 31                ', '            ', .true.,            &!431
     &  'tracer 32                ', '            ', .true.,            &!432
     &  'tracer 33                ', '            ', .true.,            &!433
     &  'tracer 34                ', '            ', .true.,            &!434
     &  'tracer 35                ', '            ', .true.,            &!435
     &  'tracer 36                ', '            ', .true.,            &!436
     &  'tracer 37                ', '            ', .true.,            &!437
     &  'tracer 38                ', '            ', .true.,            &!438
     &  'tracer 39                ', '            ', .true.,            &!439
     &  'tracer 40                ', '            ', .true.,            &!440
     &  'tracer 41                ', '            ', .true.,            &!441
     &  'tracer 42                ', '            ', .true.,            &!442
     &  'tracer 43                ', '            ', .true.,            &!443
     &  'tracer 44                ', '            ', .true.,            &!444
     &  'tracer 45                ', '            ', .true.,            &!445
     &  'tracer 46                ', '            ', .true.,            &!446
     &  'tracer 47                ', '            ', .true.,            &!447
     &  'tracer 48                ', '            ', .true.,            &!448
     &  'tracer 49                ', '            ', .true.,            &!449
     &  'tracer 50                ', '            ', .true.,            &!450
     &  'tracer 51                ', '            ', .true.,            &!451
     &  'tracer 52                ', '            ', .true.,            &!452
     &  'tracer 53                ', '            ', .true.,            &!453
     &  'tracer 54                ', '            ', .true.,            &!454
     &  'tracer 55                ', '            ', .true.,            &!455
     &  'tracer 56                ', '            ', .true.,            &!456
     &  'tracer 57                ', '            ', .true.,            &!457
     &  'tracer 58                ', '            ', .true.,            &!458
     &  'tracer 59                ', '            ', .true.,            &!459
     &  'tracer 60                ', '            ', .true.,            &!460
     &  'tracer 61                ', '            ', .true.,            &!461
     &  'tracer 62                ', '            ', .true.,            &!462
     &  'tracer 63                ', '            ', .true.,            &!463
     &  'tracer 64                ', '            ', .true.,            &!464
     &  'tracer 65                ', '            ', .true.,            &!465
     &  'tracer 66                ', '            ', .true.,            &!466
     &  'tracer 67                ', '            ', .true.,            &!467
     &  'tracer 68                ', '            ', .true.,            &!468
     &  'tracer 69                ', '            ', .true.,            &!469
     &  'tracer 70                ', '            ', .true.,            &!470
     &  'tracer 71                ', '            ', .true.,            &!471
     &  'tracer 72                ', '            ', .true.,            &!472
     &  'tracer 73                ', '            ', .true.,            &!473
     &  'tracer 74                ', '            ', .true.,            &!474
     &  'tracer 75                ', '            ', .true.,            &!475
     &  'tracer 76                ', '            ', .true.,            &!476
     &  'tracer 77                ', '            ', .true.,            &!477
     &  'tracer 78                ', '            ', .true.,            &!478
     &  'tracer 79                ', '            ', .true.,            &!479
     &  'tracer 80                ', '            ', .true.,            &!480
     &  'tracer 81                ', '            ', .true.,            &!481
     &  'tracer 82                ', '            ', .true.,            &!482
     &  'tracer 83                ', '            ', .true.,            &!483
     &  'tracer 84                ', '            ', .true.,            &!484
     &  'tracer 85                ', '            ', .true.,            &!485
     &  'tracer 86                ', '            ', .true.,            &!486
     &  'tracer 87                ', '            ', .true.,            &!487
     &  'tracer 88                ', '            ', .true.,            &!488
     &  'tracer 89                ', '            ', .true.,            &!489
     &  'tracer 90                ', '            ', .true.,            &!490
     &  'tracer 91                ', '            ', .true.,            &!491
     &  'tracer 92                ', '            ', .true.,            &!492
     &  'tracer 93                ', '            ', .true.,            &!493
     &  'tracer 94                ', '            ', .true.,            &!494
     &  'tracer 95                ', '            ', .true.,            &!495
     &  'tracer 96                ', '            ', .true.,            &!496
     &  'tracer 97                ', '            ', .true.,            &!497
     &  'tracer 98                ', '            ', .true.,            &!498
     &  'tracer 99                ', '            ', .true.,            &!499
     &  'tracer 100               ', '            ', .true./             !500

