F_V = [\
  'EGF',
  'HRG',
  'RsD',
  'RsT',
  'Kin',
  'pKin',
  'MEK',
  'A1',
  'A1_act',
  'A2',
  'A2_act',
  'A3',
  'A3_act',
  #
  'ppMEKc',
  'CREBn',
  'pCREBn',
  'ERKc',
  'ERKn',
  'pERKc',
  'pERKn',
  'ppERKc',
  'ppERKn',
  'Elk1n',
  'pElk1n',
  'cFOSc',
  'cFOSn',
  'pcFOSc',
  'pcFOSn',
  'DUSPc',
  'DUSPn',
  'pDUSPc',
  'pDUSPn',
  'DUSPn_ERKn',
  'DUSPn_pERKn',
  'DUSPn_ppERKn',
  'pDUSPn_ERKn',
  'pDUSPn_pERKn',
  'pDUSPn_ppERKn',
  'RSKc',
  'pRSKc',
  'pRSKn',
  'PrecfosmRNAn',
  'PreduspmRNAn',
  'cfosmRNAc',
  'duspmRNAc',
  #
  'len_f_vars'\
]

for i,name in enumerate(F_V):
  exec('%s=%d'%(name,i),globals())