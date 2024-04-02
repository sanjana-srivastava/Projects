function y = CDSurrogate(u)

DV =   [0.313043478	0.847826087	0.079710145	0.602173913	4.739130435	0.227536232	0.005797101
        0.473913043	0.637681159	0.442028986	0.55	4.260869565	0.21884058	0.104347826
        0.460869565	0.507246377	0.210144928	0.567391304	3.739130435	0.198550725	0.188405797
        0.47826087	0.543478261	0.405797101	0.667391304	4.347826087	0.137681159	0.084057971
        0.517391304	0.934782609	0.326086957	0.654347826	3.956521739	0.111594203	0.133333333
        0.360869565	0.768115942	0.463768116	0.67173913	4.130434783	0.105797101	0.046376812
        0.434782609	0.652173913	0.050724638	0.682608696	3.47826087	0.114492754	0.037681159
        0.5	0.927536232	0.275362319	0.626086957	3.52173913	0.285507246	0.150724638
        0.404347826	0.963768116	0.086956522	0.593478261	4.173913043	0.195652174	0.124637681
        0.495652174	0.608695652	0.31884058	0.636956522	5.260869565	0.172463768	0.168115942
        0.52173913	0.869565217	0.34057971	0.589130435	5.391304348	0.108695652	0.194202899
        0.565217391	0.797101449	0.369565217	0.576086957	4.043478261	0.282608696	0.023188406
        0.308695652	0.572463768	0.434782609	0.639130435	4.913043478	0.247826087	0.026086957
        0.426086957	0.724637681	0.384057971	0.586956522	3.043478261	0.131884058	0.095652174
        0.591304348	0.615942029	0.471014493	0.604347826	3.695652174	0.169565217	0.055072464
        0.556521739	0.52173913	0.065217391	0.584782609	3.652173913	0.117391304	0.072463768
        0.417391304	0.985507246	0.376811594	0.6	5.565217391	0.189855072	0.063768116
        0.330434783	0.891304348	0.289855072	0.676086957	5	0.288405797	0.049275362
        0.552173913	0.992753623	0.123188406	0.573913043	5.043478261	0.184057971	0.069565217
        0.395652174	0.876811594	0.297101449	0.697826087	5.608695652	0.166666667	0.12173913
        0.339130435	0.913043478	0.224637681	0.558695652	3.173913043	0.160869565	0.139130435
        0.513043478	0.586956522	0.246376812	0.610869565	3.347826087	0.128985507	0.147826087
        0.42173913	0.550724638	0.188405797	0.608695652	4.956521739	0.120289855	0.011594203
        0.443478261	0.688405797	0.173913043	0.57826087	3.304347826	0.265217391	0.101449275
        0.32173913	0.956521739	0.362318841	0.689130435	3.217391304	0.143478261	0.127536232
        0.573913043	0.630434783	0.18115942	0.569565217	5.695652174	0.155072464	0.092753623
        0.539130435	0.644927536	0.02173913	0.580434783	4.608695652	0.18115942	0.191304348
        0.369565217	0.862318841	0.427536232	0.565217391	5.130434783	0.25942029	0.04057971
        0.491304348	0.898550725	0.007246377	0.554347826	4.652173913	0.102898551	0.110144928
        0.365217391	0.804347826	0.449275362	0.680434783	4.826086957	0.250724638	0.173913043
        0.456521739	0.65942029	0.130434783	0.669565217	5.47826087	0.14057971	0.107246377
        0.486956522	0.782608696	0.239130435	0.597826087	4.434782609	0.163768116	0.086956522
        0.304347826	0.528985507	0.144927536	0.595652174	4.695652174	0.271014493	0.066666667
        0.430434783	0.84057971	0.398550725	0.693478261	3.608695652	0.294202899	0.002898551
        0.317391304	0.666666667	0.304347826	0.619565217	5.913043478	0.262318841	0.113043478
        0.560869565	0.731884058	0.195652174	0.7	5.434782609	0.279710145	0.171014493
        0.391304348	0.905797101	0.47826087	0.630434783	3.782608696	0.224637681	0.115942029
        0.413043478	0.579710145	0.347826087	0.615217391	3.913043478	0.242028986	0.142028986
        0.547826087	0.855072464	0.036231884	0.684782609	5.086956522	0.152173913	0.031884058
        0.530434783	0.594202899	0.101449275	0.591304348	4.086956522	0.230434783	0.008695652
        0.352173913	0.673913043	0.137681159	0.67826087	3.869565217	0.291304348	0.153623188
        0.595652174	0.514492754	0.15942029	0.652173913	4.217391304	0.204347826	0.130434783
        0.534782609	0.702898551	0.413043478	0.695652174	3.826086957	0.253623188	0.144927536
        0.373913043	0.536231884	0	0.686956522	4.869565217	0.244927536	0.052173913
        0.465217391	0.760869565	0.231884058	0.552173913	5.782608696	0.210144928	0.034782609
        0.452173913	0.884057971	0.253623188	0.613043478	5.52173913	0.297101449	0
        0.543478261	0.942028986	0.260869565	0.556521739	5.826086957	0.215942029	0.176811594
        0.447826087	0.623188406	0.115942029	0.643478261	4.391304348	0.1	0.2
        0.4	0.565217391	0.057971014	0.563043478	5.304347826	0.201449275	0.07826087
        0.408695652	0.746376812	0.485507246	0.663043478	5.956521739	0.256521739	0.028985507
        0.57826087	0.920289855	0.420289855	0.632608696	5.652173913	0.17826087	0.11884058
        0.347826087	0.601449275	0.5	0.658695652	3	0.123188406	0.156521739
        0.482608696	0.97826087	0.014492754	0.634782609	3.434782609	0.134782609	0.098550725
        0.326086957	0.753623188	0.202898551	0.62826087	3.391304348	0.22173913	0.185507246
        0.569565217	0.557971014	0.333333333	0.560869565	6	0.3	0.08115942
        0.37826087	1	0.072463768	0.660869565	5.347826087	0.192753623	0.017391304
        0.386956522	0.775362319	0.311594203	0.57173913	4.782608696	0.207246377	0.179710145
        0.582608696	0.68115942	0.282608696	0.691304348	4.52173913	0.236231884	0.014492754
        0.504347826	0.739130435	0.028985507	0.665217391	3.260869565	0.186956522	0.165217391
        0.508695652	0.710144928	0.043478261	0.65	4.47826087	0.273913043	0.060869565
        0.586956522	0.833333333	0.492753623	0.656521739	4.304347826	0.233333333	0.075362319
        0.334782609	0.949275362	0.217391304	0.62173913	3.086956522	0.268115942	0.057971014
        0.6	0.789855072	0.355072464	0.606521739	3.565217391	0.239130435	0.197101449
        0.526086957	0.695652174	0.456521739	0.582608696	4	0.149275362	0.182608696
        0.343478261	0.826086957	0.166666667	0.623913043	4.565217391	0.126086957	0.136231884
        0.469565217	0.811594203	0.391304348	0.645652174	5.173913043	0.146376812	0.020289855
        0.356521739	0.717391304	0.094202899	0.617391304	5.869565217	0.157971014	0.043478261
        0.382608696	0.81884058	0.152173913	0.673913043	3.130434783	0.175362319	0.089855072
        0.3	0.971014493	0.268115942	0.641304348	5.739130435	0.213043478	0.162318841
        0.439130435	0.5	0.108695652	0.647826087	5.217391304	0.276811594	0.15942029];
   
freal = [0.086106749
        0.081866643
        0.075177763
        0.085320285
        0.081978395
        0.08175073
        0.070220314
        0.065044577
        0.081645242
        0.094780999
        0.10412078
        0.073352899
        0.084643296
        0.066817828
        0.075278405
        0.077677113
        0.09974248
        0.081254021
        0.09496458
        0.098524736
        0.069404129
        0.07244752
        0.09673492
        0.065231791
        0.066801535
        0.106643916
        0.087522375
        0.08867436
        0.098117158
        0.081112967
        0.101490719
        0.087615505
        0.081474573
        0.064717573
        0.094517778
        0.085128812
        0.07263806
        0.073213362
        0.094375887
        0.077948826
        0.067814909
        0.079264021
        0.06931086
        0.083463411
        0.102231456
        0.088153783
        0.099537792
        0.089574716
        0.097085924
        0.094689492
        0.099852031
        0.065280588
        0.072043911
        0.067105666
        0.09430607
        0.095100338
        0.087749623
        0.079902977
        0.066590781
        0.077462409
        0.078046312
        0.060872155
        0.068373372
        0.08218351
        0.091729549
        0.097333046
        0.107205032
        0.064052626
        0.096077412
        0.084777978];


srgtOPT  = srgtsKRGSetOptions(DV, freal);       % Function srgtsKRGSetOptions creates the SURROGATES Toolbox option structure for kriging models
srgtOPT.KRG_RegressionModel  = @dace_regpoly2;
srgtOPT.KRG_CorrelationModel = @dace_corrgauss;
srgtOPT.KRG_Theta0           = 1e-6;
srgtOPT.KRG_LowerBound       = 1e-3;
srgtOPT.KRG_UpperBound       = 20;
srgtSRGT = srgtsKRGFit(srgtOPT);                % Function srgtsKRGFit fits the specified kriging model using the DACE toolbox of Lophaven et al. (2002).
y = srgtsKRGEvaluate(u, srgtSRGT);
end