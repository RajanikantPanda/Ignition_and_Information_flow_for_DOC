# Ignition_and_Information_flow_for_DOC
%%% For fMRI_Intrensic_Ignition_Measure
%%% Developed By: Gustavo Deco
%%% Adopted By: Rajanikant Panda 
%%% Date of Modification: 1st September 2019 
%%% Supervised: Gustavo Deco (UPF), Steven Laureys(Uliege), Jitka Annen(Uliege)  and Ane A López-González (UPF)

%%% For fMRI_Information_flow_Measure
%%% Developed By: Matthieu Gilson and Gorka Zamora-López (UPF)
%%% Adopted By: Rajanikant Panda 
%%% Date of Modification: 1st September 2019 
%%% Supervised: Matthieu Gilson and Gorka Zamora-López (UPF), Gustavo Deco (UPF), Steven Laureys(Uliege) and Jitka Annen(Uliege)

The above code compute: 
1. Intrinsic Ignition which quantify the alteration of whole-brain integration generated by naturally occurring intrinsic events
2. Tau, which measure the ‘memory depth’ of the signal. 
3. Computational model-based approach which measure the propagation of perturbations by causal spatiotemporal properties of time series given by the network effective connectivity (EC) 
   and the dynamic communicability, i.e. time-dependent graph-like descriptors. From nodal communicability, we capture the capacity of the brain region to receive and broadcast information flow after perturbations.


TO estimation of effective connectivity you need to instal the pyMOU Python package (https://github.com/mb-BCA/pyMOU) and for the study of responses to exogenous perturbations (i.e., dynamic communicability), you need to instal NetDynFlow package in Python (https://github.com/mb-BCA/NetDynFlow) developed by Matthieu Gilson and Gorka Zamora-López (UPF).
