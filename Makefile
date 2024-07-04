.PHONY: all clean

AP_MODEL = Grandi TTP Tomek ORd-Lei

# Methods
# Tuning of IKr in modified AP models
IKr_tuning:
	python3 scripts/tune_IKr.py; \
	python3 scripts/fig_scripts/IKr_tuning.py; \

# Background figures
# Current contribution of AP models
# Require conductance scaling factor to be calculated
APmodels:
	python3 scripts/fig_scripts/base_APmodels.py; \

# State occupancy of IKr models
IKrmodels:
	python3 scripts/base_hERG.py; \
	python3 scripts/fig_scripts/base_hERGmodels.py; \

# Results
# APD90 and qNet comparison between dofetilide and verapamil for all AP models
dof_vs_ver:
	python3 scripts/compare_binding_kinetics_all.py; \
	python3 scripts/fig_scripts/compile_AP_compare.py; \

# RMSDs between the SD model and the CS model for all AP models
param_space_exploration:
	for model in $(AP_MODEL); do \
		python3 scripts/create_param_space.py $$model; \
		python3 scripts/SA_param_space.py $$model; \
		python3 scripts/combine_results.py $$model; \
	done
	python3 scripts/fig_scripts/SA_3D_all.py; \
	python3 scripts/fig_scripts/SA_3D.py ORd-Lei; \

# SA.sh

# Supplementary
# RMSD distribution with varying Hill coefficient
# To define similarity range
SA_drugs_n:
	for model in $(AP_MODEL); do \
		python3 script/SA_param_interest.py $$model; \
	done
	python3 script/fig_scripts/SA_N_diff.py; \

# Optimisation of parameters for Lei-SD model
param_opt:
	python3 scripts/side_scripts/parameter_fit.py; \
	python3 scripts/side_scripts/compile_param_fit.py; \
	python3 scripts/side_scripts/param_fit_fig.py; \

# Extra
# AP comparison between the SD model and the CS model for each AP models
AP_compare:
	for drug in dofetilide verapamil; do \
		for model in $(AP_MODEL); do \
			python3 scripts/fig_scripts/APD_compare.py $$model $$drug; \
		done \
	done
