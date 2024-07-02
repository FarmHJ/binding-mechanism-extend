.PHONY: all clean

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
# Biomarker comparison between dofetilide and verapamil for all AP models
dof_vs_ver:
	for model in $(AP_MODEL); do \
		python3 scripts/compare_binding_kinetics.py $$model dofetilide --ikr_tuning AP_duration; \
		python3 scripts/compare_binding_kinetics.py $$model verapamil --ikr_tuning AP_duration; \
	done
	python3 scripts/fig_scripts/compile_AP_compare.py; \

param_space_exploration:
	for model in $(AP_MODEL); do \
		python3 scripts/SA_param_space.py $$model --ikr_tuning AP_duration; \
		python3 scripts/fig_scripts/SA_3D.py $$model; \
	done

# Supplementary
# Optimisation of parameters for Lei-SD model
param_opt:
	python3 scripts/side_scripts/parameter_fit.py; \
	python3 scripts/side_scripts/compile_param_fit.py; \
	python3 scripts/side_scripts/param_fit_fig.py; \

# RMSD distribution with varying Hill coefficient
SA_Hill:
	for model in TTP Tomek ORd-Lei; do \
		python3 scripts/SA_param_interest.py $$model --ikr_tuning AP_duration; \
		python3 scripts/SA_param_interest.py $$model --mode parameter_SA --ikr_tuning AP_duration; \
		python3 scripts/fig_scripts/SA_N_diff.py $$model; \
	done

# Biomarker comparison between dofetilide and verapamil for each AP models
dof_vs_ver_AP:
	for drug in dofetilide verapamil; do \
		for model in $(AP_MODEL); do \
			python3 scripts/compare_binding_kinetics.py $$model $$drug --ikr_tuning AP_duration -m AP; \
			python3 scripts/fig_scripts/APD_compare.py $$model $$drug --ikr_tuning AP_duration; \
		done \
	done
