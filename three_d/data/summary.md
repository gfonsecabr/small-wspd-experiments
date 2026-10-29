# Experiments in R^3: numbers for the text (analysis/exp3d.py; seed 1 of the generators and the scans)

instances per n: 200: 10, 300: 10, 400: 10, 500: 10, 600: 10, 700: 10, 800: 10, 900: 10, 1000: 10, 2000: 10, 3000: 10, 4000: 10, 5000: 10, 6000: 10, 7000: 10, 8000: 10, 9000: 10, 10000: 10

## size / reference (mean, max) over all instances
n=200 (ref LB on 0): CK box 1.0749 / 1.3654; CK 1.0586 / 1.2239; Sweep 1.0084 / 1.0689; CK + reshape 1.0000 / 1.0000; Sweep + reshape 1.0000 / 1.0003
n=300 (ref LB on 0): CK box 1.1160 / 1.5818; CK 1.0910 / 1.3724; Sweep 1.0103 / 1.0753; CK + reshape 1.0000 / 1.0001; Sweep + reshape 1.0000 / 1.0003
n=400 (ref LB on 0): CK box 1.1207 / 1.4835; CK 1.0923 / 1.2590; Sweep 1.0117 / 1.0804; CK + reshape 1.0000 / 1.0000; Sweep + reshape 1.0000 / 1.0001
n=500 (ref LB on 0): CK box 1.1520 / 1.6258; CK 1.1178 / 1.3607; Sweep 1.0141 / 1.0934; CK + reshape 1.0000 / 1.0003; Sweep + reshape 1.0000 / 1.0001
n=600 (ref LB on 0): CK box 1.1645 / 1.6595; CK 1.1289 / 1.3877; Sweep 1.0170 / 1.1054; CK + reshape 1.0000 / 1.0002; Sweep + reshape 1.0000 / 1.0001
n=700 (ref LB on 0): CK box 1.1786 / 1.6716; CK 1.1414 / 1.3932; Sweep 1.0179 / 1.1003; CK + reshape 1.0000 / 1.0001; Sweep + reshape 1.0000 / 1.0002
n=800 (ref LB on 0): CK box 1.1939 / 1.7180; CK 1.1534 / 1.4198; Sweep 1.0205 / 1.1108; CK + reshape 1.0000 / 1.0001; Sweep + reshape 1.0000 / 1.0001
n=900 (ref LB on 0): CK box 1.2050 / 1.6970; CK 1.1653 / 1.4272; Sweep 1.0222 / 1.1152; CK + reshape 1.0000 / 1.0002; Sweep + reshape 1.0000 / 1.0002
n=1000 (ref LB on 0): CK box 1.2201 / 1.7604; CK 1.1775 / 1.4753; Sweep 1.0240 / 1.1170; CK + reshape 1.0000 / 1.0001; Sweep + reshape 1.0000 / 1.0002
n=2000 (ref LB on 10): CK box 1.3038 / 1.8640; CK 1.2455 / 1.5417; Sweep 1.0421 / 1.1407; CK + reshape 1.0042 / 1.0156; Sweep + reshape 1.0042 / 1.0156
n=3000 (ref LB on 10): CK box 1.3672 / 1.9092; CK 1.2953 / 1.5572; Sweep 1.0560 / 1.1480; CK + reshape 1.0061 / 1.0207; Sweep + reshape 1.0061 / 1.0206
n=4000 (ref LB on 10): CK box 1.4196 / 1.9414; CK 1.3383 / 1.5873; Sweep 1.0688 / 1.1570; CK + reshape 1.0079 / 1.0234; Sweep + reshape 1.0079 / 1.0234
n=5000 (ref LB on 10): CK box 1.4603 / 1.9696; CK 1.3694 / 1.6067; Sweep 1.0805 / 1.1610; CK + reshape 1.0098 / 1.0244; Sweep + reshape 1.0098 / 1.0244
n=6000 (ref LB on 10): CK box 1.5018 / 1.9844; CK 1.4011 / 1.6080; Sweep 1.0921 / 1.1654; CK + reshape 1.0118 / 1.0253; Sweep + reshape 1.0118 / 1.0254
n=7000 (ref LB on 10): CK box 1.5314 / 1.9924; CK 1.4246 / 1.6233; Sweep 1.1021 / 1.1645; CK + reshape 1.0140 / 1.0258; Sweep + reshape 1.0140 / 1.0259
n=8000 (ref LB on 10): CK box 1.5638 / 2.0175; CK 1.4494 / 1.6428; Sweep 1.1121 / 1.1887; CK + reshape 1.0160 / 1.0299; Sweep + reshape 1.0161 / 1.0300
n=9000 (ref LB on 10): CK box 1.5879 / 2.0204; CK 1.4686 / 1.6586; Sweep 1.1217 / 1.2114; CK + reshape 1.0183 / 1.0366; Sweep + reshape 1.0184 / 1.0367
n=10000 (ref LB on 10): CK box 1.6115 / 2.0341; CK 1.4864 / 1.6840; Sweep 1.1311 / 1.2327; CK + reshape 1.0209 / 1.0432; Sweep + reshape 1.0209 / 1.0434

## per input at n = 10000 (size / reference)
armadillo: CK box 1.832; CK 1.668; Sweep 1.202; CK + reshape 1.035; Sweep + reshape 1.035; pairs/dumbbell 9.01; dumbbells/point (sweep+reshape) 574.7
buddha: CK box 1.757; CK 1.606; Sweep 1.159; CK + reshape 1.023; Sweep + reshape 1.023; pairs/dumbbell 7.60; dumbbells/point (sweep+reshape) 672.4
bunny: CK box 1.807; CK 1.667; Sweep 1.233; CK + reshape 1.043; Sweep + reshape 1.043; pairs/dumbbell 8.25; dumbbells/point (sweep+reshape) 632.0
cluster: CK box 2.034; CK 1.656; Sweep 1.169; CK + reshape 1.015; Sweep + reshape 1.015; pairs/dumbbell 14.68; dumbbells/point (sweep+reshape) 345.6
dragon: CK box 1.779; CK 1.621; Sweep 1.166; CK + reshape 1.026; Sweep + reshape 1.026; pairs/dumbbell 8.47; dumbbells/point (sweep+reshape) 606.2
gauss: CK box 1.296; CK 1.267; Sweep 1.028; CK + reshape 1.003; Sweep + reshape 1.003; pairs/dumbbell 1.90; dumbbells/point (sweep+reshape) 2645.9
grid: CK box 1.052; CK 1.051; Sweep 1.019; CK + reshape 1.007; Sweep + reshape 1.007; pairs/dumbbell 1.16; dumbbells/point (sweep+reshape) 4347.5
shell: CK box 1.871; CK 1.684; Sweep 1.182; CK + reshape 1.026; Sweep + reshape 1.026; pairs/dumbbell 5.96; dumbbells/point (sweep+reshape) 861.1
spiral: CK box 1.319; CK 1.317; Sweep 1.120; CK + reshape 1.028; Sweep + reshape 1.028; pairs/dumbbell 243.83; dumbbells/point (sweep+reshape) 21.1
uniform: CK box 1.368; CK 1.326; Sweep 1.033; CK + reshape 1.003; Sweep + reshape 1.003; pairs/dumbbell 2.23; dumbbells/point (sweep+reshape) 2252.1

## smallest construction
n=200: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=300: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=400: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=500: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=600: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=700: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=800: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=900: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=1000: sweep <= min(CK, CK box) on 10 of 10, < on 9; CK < sweep on 0 ()
n=2000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=3000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=4000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=5000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=6000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=7000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=8000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=9000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()
n=10000: sweep <= min(CK, CK box) on 10 of 10, < on 10; CK < sweep on 0 ()

## reshape: sweep start / CK start (final sizes); against OPT
all: max |1 - sweep/ck| 0.00026
lsck on proven OPT (90): mean 1.00002 max 1.00032, = OPT on 72
lssweep on proven OPT (90): mean 1.00002 max 1.00027, = OPT on 69

## lower bounds
n=200 proven=10: LB0/OPT mean 0.9880 min 0.9198; LB1/OPT 0.9963; LB/OPT mean 0.99928 min 0.99476, = OPT on 6; alpha/OPT mean 0.99953 min 0.99581 (N=10), LB = alpha on 8
n=300 proven=10: LB0/OPT mean 0.9859 min 0.9159; LB1/OPT 0.9951; LB/OPT mean 0.99879 min 0.99330, = OPT on 3; alpha/OPT mean 0.99935 min 0.99491 (N=10), LB = alpha on 8
n=400 proven=10: LB0/OPT mean 0.9831 min 0.9089; LB1/OPT 0.9958; LB/OPT mean 0.99868 min 0.99186, = OPT on 3; alpha/OPT mean 0.99906 min 0.99218 (N=10), LB = alpha on 8
n=500 proven=10: LB0/OPT mean 0.9788 min 0.8883; LB1/OPT 0.9943; LB/OPT mean 0.99790 min 0.98676, = OPT on 2; alpha/OPT mean 0.99857 min 0.98764 (N=10), LB = alpha on 7
n=600 proven=10: LB0/OPT mean 0.9772 min 0.8910; LB1/OPT 0.9937; LB/OPT mean 0.99789 min 0.98930, = OPT on 2; alpha/OPT mean nan min nan (N=0), LB = alpha on 0
n=700 proven=10: LB0/OPT mean 0.9754 min 0.8895; LB1/OPT 0.9936; LB/OPT mean 0.99783 min 0.98909, = OPT on 1; alpha/OPT mean nan min nan (N=0), LB = alpha on 0
n=800 proven=10: LB0/OPT mean 0.9726 min 0.8824; LB1/OPT 0.9931; LB/OPT mean 0.99771 min 0.98955, = OPT on 1; alpha/OPT mean nan min nan (N=0), LB = alpha on 0
n=900 proven=10: LB0/OPT mean 0.9691 min 0.8662; LB1/OPT 0.9911; LB/OPT mean 0.99741 min 0.98753, = OPT on 1; alpha/OPT mean nan min nan (N=0), LB = alpha on 0
n=1000 proven=10: LB0/OPT mean 0.9687 min 0.8795; LB1/OPT 0.9915; LB/OPT mean 0.99766 min 0.99073, = OPT on 1; alpha/OPT mean nan min nan (N=0), LB = alpha on 0
alpha instances 40: of the total gap alpha - LB0, swaps close 0.853 and plateau moves 0.137
n=200: LB/LB0 1.0120 (max 1.0815); sweep+reshape / LB 1.0007 (max 1.0055)
n=300: LB/LB0 1.0137 (max 1.0845); sweep+reshape / LB 1.0012 (max 1.0070)
n=400: LB/LB0 1.0165 (max 1.0913); sweep+reshape / LB 1.0013 (max 1.0083)
n=500: LB/LB0 1.0205 (max 1.1109); sweep+reshape / LB 1.0021 (max 1.0135)
n=600: LB/LB0 1.0222 (max 1.1103); sweep+reshape / LB 1.0021 (max 1.0109)
n=700: LB/LB0 1.0239 (max 1.1119); sweep+reshape / LB 1.0022 (max 1.0112)
n=800: LB/LB0 1.0269 (max 1.1215); sweep+reshape / LB 1.0023 (max 1.0106)
n=900: LB/LB0 1.0306 (max 1.1400); sweep+reshape / LB 1.0026 (max 1.0128)
n=1000: LB/LB0 1.0311 (max 1.1266); sweep+reshape / LB 1.0024 (max 1.0094)
n=2000: LB/LB0 1.0461 (max 1.1475); sweep+reshape / LB 1.0042 (max 1.0156)
n=3000: LB/LB0 1.0569 (max 1.1604); sweep+reshape / LB 1.0061 (max 1.0206)
n=4000: LB/LB0 1.0647 (max 1.1636); sweep+reshape / LB 1.0079 (max 1.0234)
n=5000: LB/LB0 1.0703 (max 1.1620); sweep+reshape / LB 1.0098 (max 1.0244)
n=6000: LB/LB0 1.0757 (max 1.1665); sweep+reshape / LB 1.0118 (max 1.0254)
n=7000: LB/LB0 1.0803 (max 1.1678); sweep+reshape / LB 1.0140 (max 1.0259)
n=8000: LB/LB0 1.0846 (max 1.1710); sweep+reshape / LB 1.0161 (max 1.0300)
n=9000: LB/LB0 1.0881 (max 1.1727); sweep+reshape / LB 1.0184 (max 1.0367)
n=10000: LB/LB0 1.0912 (max 1.1752); sweep+reshape / LB 1.0209 (max 1.0434)
n=200: time sweep 0.02 s, LB1 after the sweep 0.00 s, LB after the sweep 0.02 s
n=300: time sweep 0.05 s, LB1 after the sweep 0.01 s, LB after the sweep 0.04 s
n=400: time sweep 0.08 s, LB1 after the sweep 0.01 s, LB after the sweep 0.08 s
n=500: time sweep 0.13 s, LB1 after the sweep 0.02 s, LB after the sweep 0.15 s
n=600: time sweep 0.18 s, LB1 after the sweep 0.03 s, LB after the sweep 0.21 s
n=700: time sweep 0.24 s, LB1 after the sweep 0.04 s, LB after the sweep 0.30 s
n=800: time sweep 0.32 s, LB1 after the sweep 0.06 s, LB after the sweep 0.41 s
n=900: time sweep 0.41 s, LB1 after the sweep 0.08 s, LB after the sweep 0.54 s
n=1000: time sweep 0.50 s, LB1 after the sweep 0.10 s, LB after the sweep 0.68 s
n=2000: time sweep 1.98 s, LB1 after the sweep 0.45 s, LB after the sweep 3.37 s
n=3000: time sweep 4.27 s, LB1 after the sweep 1.14 s, LB after the sweep 8.65 s
n=4000: time sweep 7.12 s, LB1 after the sweep 2.16 s, LB after the sweep 17.75 s
n=5000: time sweep 10.85 s, LB1 after the sweep 3.72 s, LB after the sweep 32.79 s
n=6000: time sweep 15.80 s, LB1 after the sweep 6.62 s, LB after the sweep 69.52 s
n=7000: time sweep 20.79 s, LB1 after the sweep 9.88 s, LB after the sweep 94.33 s
n=8000: time sweep 26.70 s, LB1 after the sweep 14.28 s, LB after the sweep 131.07 s
n=9000: time sweep 34.04 s, LB1 after the sweep 18.11 s, LB after the sweep 174.65 s
n=10000: time sweep 41.63 s, LB1 after the sweep 26.13 s, LB after the sweep 233.13 s

## OPT and alpha (CPLEX)
n=200: proven 10, time median 0.03 s mean 0.2 s max 1 s
n=300: proven 10, time median 0.13 s mean 0.9 s max 6 s
n=400: proven 10, time median 0.19 s mean 3.0 s max 24 s
n=500: proven 10, time median 0.37 s mean 8.6 s max 75 s
n=600: proven 10, time median 0.63 s mean 18.7 s max 171 s
n=700: proven 10, time median 1.05 s mean 50.3 s max 477 s
n=800: proven 10, time median 1.63 s mean 104.6 s max 1015 s
n=900: proven 10, time median 2.41 s mean 182.0 s max 1776 s
n=1000: proven 10, time median 3.18 s mean 333.7 s max 3277 s
alpha: 40 runs, optimal 40, time max 88 s

## running times (mean s) and peak memory (max GB)
n=200: CK box 0.00; CK 0.00; Sweep 0.02
n=300: CK box 0.00; CK 0.00; Sweep 0.05
n=400: CK box 0.00; CK 0.00; Sweep 0.08
n=500: CK box 0.01; CK 0.01; Sweep 0.13
n=600: CK box 0.01; CK 0.01; Sweep 0.19
n=700: CK box 0.01; CK 0.01; Sweep 0.25
n=800: CK box 0.02; CK 0.02; Sweep 0.34
n=900: CK box 0.02; CK 0.02; Sweep 0.43
n=1000: CK box 0.02; CK 0.03; Sweep 0.53
n=2000: CK box 0.09; CK 0.10; Sweep 2.11
n=3000: CK box 0.17; CK 0.19; Sweep 4.55
n=4000: CK box 0.29; CK 0.31; Sweep 7.67
n=5000: CK box 0.43; CK 0.46; Sweep 11.61
n=6000: CK box 0.64; CK 0.68; Sweep 17.11
n=7000: CK box 0.79; CK 0.83; Sweep 22.54
n=8000: CK box 0.95; CK 1.03; Sweep 28.80
n=9000: CK box 1.21; CK 1.30; Sweep 36.91
n=10000: CK box 1.44; CK 1.50; Sweep 45.03
memory n=10000 ck: max 5.01 GB
memory n=10000 ckexact: max 5.01 GB
memory n=10000 lb_swap3: max 0.88 GB
memory n=10000 lsck: max 8.37 GB
memory n=10000 lssweep: max 8.25 GB
memory n=10000 sweep: max 8.42 GB

## local search on the samples (fig:ls3d): mean size / reference at the end of the budget
n=200 lns from ck (10): 1.00080
n=200 reshape from ck (10): 1.00000
n=200 lns from sweep (10): 1.00049
n=200 reshape from sweep (10): 1.00003
n=1000 lns from ck (10): 1.00363
n=1000 reshape from ck (10): 1.00002
n=1000 lns from sweep (10): 1.00192
n=1000 reshape from sweep (10): 1.00003
n=10000 lns from ck (3): 1.33512
n=10000 reshape from ck (3): 1.02029
n=10000 lns from sweep (3): 1.10117
n=10000 reshape from sweep (3): 1.02035
n=200 from ck: reshape <= LNS on 10 of 10, equal on 9
n=200 from sweep: reshape <= LNS on 10 of 10, equal on 8
n=1000 from ck: reshape <= LNS on 10 of 10, equal on 2
n=1000 from sweep: reshape <= LNS on 10 of 10, equal on 2
n=10000 from ck: reshape <= LNS on 3 of 3, equal on 0
n=10000 from sweep: reshape <= LNS on 3 of 3, equal on 0

## partition WSPDs (size / reference, mean)
n=200: rs 1.0000; rs_pn 1.0025; rs_p 1.0015; rc_p 1.0017; sw_pn 1.0113; sw_p 1.0088; cke 1.0586; cke_m 1.0162; qte 1.2417; qte_m 1.0243; N=10
n=300: rs 1.0000; rs_pn 1.0030; rs_p 1.0019; rc_p 1.0021; sw_pn 1.0142; sw_p 1.0107; cke 1.0910; cke_m 1.0225; qte 1.2708; qte_m 1.0285; N=10
n=400: rs 1.0000; rs_pn 1.0037; rs_p 1.0025; rc_p 1.0026; sw_pn 1.0164; sw_p 1.0124; cke 1.0923; cke_m 1.0263; qte 1.3247; qte_m 1.0385; N=10
n=500: rs 1.0000; rs_pn 1.0050; rs_p 1.0031; rc_p 1.0034; sw_pn 1.0201; sw_p 1.0149; cke 1.1178; cke_m 1.0367; qte 1.3699; qte_m 1.0561; N=10
n=600: rs 1.0000; rs_pn 1.0056; rs_p 1.0036; rc_p 1.0030; sw_pn 1.0250; sw_p 1.0188; cke 1.1289; cke_m 1.0458; qte 1.4127; qte_m 1.0717; N=10
n=700: rs 1.0000; rs_pn 1.0060; rs_p 1.0040; rc_p 1.0039; sw_pn 1.0253; sw_p 1.0190; cke 1.1414; cke_m 1.0550; qte 1.4434; qte_m 1.0889; N=10
n=800: rs 1.0000; rs_pn 1.0064; rs_p 1.0045; rc_p 1.0045; sw_pn 1.0294; sw_p 1.0220; cke 1.1534; cke_m 1.0654; qte 1.4631; qte_m 1.1062; N=10
n=900: rs 1.0000; rs_pn 1.0071; rs_p 1.0047; rc_p 1.0048; sw_pn 1.0319; sw_p 1.0243; cke 1.1653; cke_m 1.0744; qte 1.4832; qte_m 1.1222; N=10
n=1000: rs 1.0000; rs_pn 1.0071; rs_p 1.0049; rc_p 1.0050; sw_pn 1.0345; sw_p 1.0264; cke 1.1775; cke_m 1.0853; qte 1.5122; qte_m 1.1434; N=10
n=2000: rs 1.0042; rs_pn 1.0158; rs_p 1.0128; rc_p 1.0127; sw_pn 1.0617; sw_p 1.0505; cke 1.2455; cke_m 1.1765; qte 1.7318; qte_m 1.3592; N=10
n=3000: rs 1.0061; rs_pn 1.0214; rs_p 1.0179; rc_p 1.0178; sw_pn 1.0836; sw_p 1.0705; cke 1.2953; cke_m 1.2376; qte 1.8666; qte_m 1.5336; N=10
n=4000: rs 1.0079; rs_pn 1.0268; rs_p 1.0229; rc_p 1.0226; sw_pn 1.1061; sw_p 1.0914; cke 1.3383; cke_m 1.2896; qte 1.9750; qte_m 1.6714; N=10
n=5000: rs 1.0098; rs_pn 1.0323; rs_p 1.0278; rc_p 1.0276; sw_pn 1.1275; sw_p 1.1116; cke 1.3694; cke_m 1.3263; qte 2.0671; qte_m 1.7960; N=10
n=6000: rs 1.0118; rs_pn 1.0374; rs_p 1.0326; rc_p 1.0327; sw_pn 1.1493; sw_p 1.1323; cke 1.4011; cke_m 1.3624; qte 2.1646; qte_m 1.9123; N=10
n=7000: rs 1.0140; rs_pn 1.0430; rs_p 1.0379; rc_p 1.0378; sw_pn 1.1706; sw_p 1.1524; cke 1.4246; cke_m 1.3901; qte 2.2439; qte_m 2.0128; N=10
n=8000: rs 1.0161; rs_pn 1.0486; rs_p 1.0431; rc_p 1.0429; sw_pn 1.1920; sw_p 1.1729; cke 1.4494; cke_m 1.4166; qte 2.3155; qte_m 2.1025; N=10
n=9000: rs 1.0184; rs_pn 1.0551; rs_p 1.0491; rc_p 1.0489; sw_pn 1.2136; sw_p 1.1936; cke 1.4686; cke_m 1.4380; qte 2.3809; qte_m 2.1811; N=10
n=10000: rs 1.0209; rs_pn 1.0615; rs_p 1.0551; rc_p 1.0548; sw_pn 1.2354; sw_p 1.2145; cke 1.4864; cke_m 1.4576; qte 2.4392; qte_m 2.2479; N=10
converted sweep + reshape (+ merge) larger than CK + merge on 0 of 180; than CK on 0
octree + merge < CK + merge on 5 of 180
