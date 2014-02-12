

integer,parameter :: MAXMP=35
integer,parameter :: Ind_00=1,Ind_10=2,Ind_01=3,Ind_20=4,Ind_02=5, &
                     Ind_11=6,Ind_30=7,Ind_03=8,Ind_21=9,Ind_12=10
integer,parameter :: Ind_000=1,Ind_100=2,Ind_010=3,Ind_001=4,Ind_200=5, &
                     Ind_020=6,Ind_002=7,Ind_110=8,Ind_101=9,Ind_011=10
integer,parameter :: Ind_300=11,Ind_030=12,Ind_003=13,Ind_210=14,Ind_201=15,&
                     Ind_120=16,Ind_021=17,Ind_102=18,Ind_012=19,Ind_111=20
integer,parameter :: Ind_400=21,Ind_040=22,Ind_004=23,Ind_310=24,Ind_301=25, &
                     Ind_130=26,Ind_031=27,Ind_103=28,Ind_013=29,Ind_220=30, &
                     Ind_202=31,Ind_022=32,Ind_211=33,Ind_121=34,Ind_112=35
integer,parameter :: Ind_500=36,Ind_050=37,Ind_005=38,Ind_410=39,Ind_401=40, &
                     Ind_140=41,Ind_041=42,Ind_104=43,Ind_014=44,Ind_320=45, &
                     Ind_302=46,Ind_230=47,Ind_032=48,Ind_203=49,Ind_023=50, &
                     Ind_311=51,Ind_131=52,Ind_113=53,Ind_221=54,Ind_212=55, &
                     Ind_122=56

integer,dimension(10),parameter :: deriv1 =  &
                                        (/Ind_100,Ind_200,Ind_110, &
                                          Ind_101,Ind_300,Ind_120, &
                                          Ind_102,Ind_210,Ind_201,Ind_111/)
integer,dimension(10),parameter :: deriv2 =  &
                                        (/Ind_010,Ind_110,Ind_020, &
                                          Ind_011,Ind_210,Ind_030, &
                                          Ind_012,Ind_120,Ind_111,Ind_021/)
integer,dimension(10),parameter :: deriv3 =  &
                                        (/Ind_001,Ind_101,Ind_011, &
                                          Ind_002,Ind_201,Ind_021, &
                                          Ind_003,Ind_111,Ind_102,Ind_012/)
