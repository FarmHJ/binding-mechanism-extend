.PHONY: all clean

AP_MODELS = Grandi TTP Tomek-Cl
DRUGS = dofetilide verapamil

# Get Hill coef of drugs for Lei model
Hill_Lei:
	cd scripts/; \
	for drug in $(DRUGS); do\
		python3 drug_effect_Lei.py $$drug --plot; \
	done

# Get IKr conductance scale for all AP-IKr models
tune_IKr:
	cd scripts/; \
	for ap_model in $(AP_MODELS); do \
		python3 tune_IKr.py $$ap_model; \
	done
	cd scripts/; \
	python3 tune_IKr.py ORd-Lei;

# Compare kinetics
kinetics:
	# cd scripts/; \
	# for ap_model in $(AP_MODELS); do \
	# 	python3 compare_binding_kinetics.py $$ap_model dofetilide; \
	# 	python3 compare_binding_kinetics.py $$ap_model verapamil; \
	# done
	cd scripts/; \
	python3 compare_binding_kinetics.py ORd-Lei dofetilide --ikr_tuning AP_duration; \
	python3 compare_binding_kinetics.py ORd-Lei verapamil --ikr_tuning AP_duration; \

kinetics_apd90:
	cd scripts/; \
	for ap_model in $(AP_MODELS); do \
		python3 compare_binding_kinetics.py $$ap_model dofetilide -m AP; \
		python3 compare_binding_kinetics.py $$ap_model verapamil -m AP; \
	done
	cd scripts/; \
	python3 compare_binding_kinetics.py ORd-Lei dofetilide --ikr_tuning AP_duration -m AP; \
	python3 compare_binding_kinetics.py ORd-Lei verapamil --ikr_tuning AP_duration -m AP; \

kinetics_fig:
	# cd scripts/fig_scripts/; \
	# for ap_model in $(AP_MODELS); do \
	# 	python3 APD_compare.py $$ap_model dofetilide; \
	# 	python3 APD_compare.py $$ap_model verapamil; \
	# done
	cd scripts/fig_scripts/; \
	python3 APD_compare.py ORd-Lei dofetilide --ikr_tuning AP_duration; \
	python3 APD_compare.py ORd-Lei verapamil --ikr_tuning AP_duration; \

kinetics_fig_all:
	cd scripts/fig_scripts/; \
	python3 compile_AP_compare.py

RMSE_drugs:
	# cd scripts/; \
	# for ap_model in $(AP_MODELS); do \
	# 	python3 SA_param_interest.py $$ap_model --mode only_drugs; \
	# done
	cd scripts/; \
	python3 SA_param_interest.py ORd-Lei --mode only_drugs --ikr_tuning AP_duration; \

SA_n:
	# cd scripts/; \
	# for ap_model in $(AP_MODELS); do \
	# 	python3 SA_param_interest.py $$ap_model --mode parameter_SA; \
	# done
	cd scripts/; \
	python3 SA_param_interest.py ORd-Lei --mode parameter_SA --ikr_tuning AP_duration; \

AP_schematic:
	cd scripts/; \
	python3 side_scripts/base_hERG.py; \
	python3 fig_scripts/AP_schematic.py; \

DRUG_LIST = bepridil terfenadine cisapride ranolazine quinidine sotalol chlorpromazine ondansetron diltiazem mexiletine
AP_MODEL = ORd-Li Grandi TTP Tomek
IKR_MODEL = Li Lei

# TODO: check the data required for load_Hill_eq function
# TODO: rerun everything with AP_duration tuning method

# Background figures
# Current contribution of AP models
# Check the IKr_tuning method
APmodels:
	python3 scripts/fig_scripts/base_APmodels.py; \

# State occupancy of IKr models
IKrmodels:
	python3 scripts/base_hERG.py; \
	python3 scripts/fig_scripts/base_hERGmodels.py; \

# Methods
# Tuning of IKr in modified AP models
IKr_tuning:
	for model in $(AP_MODEL); do \
		python3 scripts/tune_IKr.py $$model --method AP_duration --noplot; \
	done
	python3 scripts/fig_scripts/IKr_tuning.py; \

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
	python3 scripts/side_scripts/param_fit_fig.py; \
	python3 scripts/side_scripts/compile_param_fit.py; \

# RMSD distribution with varying Hill coefficient
SA_Hill:
	for model in $(AP_MODEL); do \
		python3 scripts/SA_param_interest.py $$model --ikr_tuning AP_duration; \
		python3 scripts/SA_param_interest.py $$model --mode parameter_SA --ikr_tuning AP_duration; \
		python3 scripts/fig_scripts/SA_N_diff.py $$model
	done

# Biomarker comparison between dofetilide and verapamil for each AP models
dof_vs_ver_AP:
	for model in $(AP_MODEL); do \
		python3 scripts/fig_scripts/APD_compare.py $$model dofetilide --ikr_tuning AP_duration; \
		python3 scripts/fig_scripts/APD_compare.py $$model verapamil --ikr_tuning AP_duration; \
	done
