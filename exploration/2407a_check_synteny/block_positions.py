# %%
import pypangraph as pp
import pandas as pd

fname = "../../results/ST131_ABC/pangraph/asm20-100-5-polished.json"
pan = pp.Pangraph.load_json(fname)

# LATJWEMOFA,222,8484,#ff2100
# WDPQHEJPPO,222,185990,#f10000
# BKGMWIBGXJ,222,155008,#db0000
# WDXIDBBCPP,222,99842,#d00000
# XXVMWZCEKI,222,125489,#cc4c4c
# HEVMPTBCLE,222,1340,#cccccc
# XISLAJEPQB,222,495997,#000000
# RVLRLPDSYQ,222,67421,#4b0055
# QQSILILDBT,222,36253,#7b008c
# ELYKQDZNBM,222,827,#860097
# KGJXRPELGX,222,1919,#3800a3
# KIVDDRSTJR,222,17704,#0000b5
# SKCSCYCISB,222,11874,#0000d5
# MVMOFPVELT,222,28685,#0038dd
# WNLYKTZGKW,222,18069,#007ddd
# QJRASUKHLX,222,3426,#0092dd
# IFRPFFEGON,222,105716,#00a0c7
# GPXHVRRZLC,222,408241,#00aaa8
# KYQOKYBCOW,222,60067,#00aa90
# YBWQKVQGZE,222,76973,#00a353
# XRXZJDDTTM,222,140178,#009a00
# EQHTXGCHGZ,222,623,#00af00
# BLEICFLMSM,222,3572,#00c700
# NEECSYVOPQ,222,5865,#00dc00
# KPBYGJHRZJ,222,293868,#00f200
# CYGJWOEQKN,222,3823,#2cff00
# TORJAESOXF,222,477776,#b0ff00
# RYYAQMEJGY,222,114706,#d8f500
# PLTCZQCVRD,222,112403,#f1e700
# ZPIAIJAZDD,222,33819,#fcd200
# KSXZRIJFEH,222,64853,#ffb100
# KKPYPKGMXA,222,407537,#ff8100

# %%
idxs = [
    ("QQSILILDBT", "NZ_CP049085.2"),
    ("KKPYPKGMXA", "NZ_CP124320.1"),
    ("GPXHVRRZLC", "NZ_CP124328.1"),
    ("KIVDDRSTJR", "NZ_CP124328.1"),
    ("QJRASUKHLX", "NZ_CP124372.1"),
    ("CYGJWOEQKN", "NZ_OX637959.1"),
    ("CYGJWOEQKN", "NZ_OX637961.1"),
    ("WDPQHEJPPO", "NZ_CP124481.1"),
    ("XRXZJDDTTM", "NZ_OW968277.1"),
]

for b, iso in idxs:
    p = pan.paths[iso]
    b_idx = list(p.block_ids).index(b)
    s, e = p.block_positions[b_idx : b_idx + 2]
    o = p.block_strands[b_idx]
    print(f"{b=}, {iso=}, {s=}, {e=}, {o=}")
# %%
