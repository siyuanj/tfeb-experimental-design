"""Shared NPG/nrc colors from ggsci; neutral text and white backgrounds."""

SOURCE = 'https://github.com/nanxstats/ggsci/blob/37b74f85e62680b9d4523b0b4c0d9bfa0403d299/R/palettes.R'
NPG = ('#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F',
       '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85')
INK = '#222222'
BLUE, TEAL, ORANGE = NPG[3], NPG[2], NPG[0]
FEATURE_COLORS = {'CMV promoter': BLUE, 'TFEB': TEAL, 'EGFP': NPG[1],
                  'KanR/NeoR': NPG[5], 'pUC origin': NPG[9], 'f1 origin': NPG[8],
                  'linker GDPPVAT': NPG[4], 'TFEB-EGFP fusion ORF': BLUE}
