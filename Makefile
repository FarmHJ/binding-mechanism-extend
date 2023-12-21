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

