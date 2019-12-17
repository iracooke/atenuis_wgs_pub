for model in 'no_mig' 'sym_mig' 'asym_mig' 'priorsize_asym_mig' 'asym_mig_size' 'isolation_asym_mig';do 
	cat optim/${model}_[0-9].${model}.optimized.txt | grep -v 'Model' | awk -f best_model.awk 
done

