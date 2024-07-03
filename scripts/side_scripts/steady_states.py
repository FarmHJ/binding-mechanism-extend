import myokit
import os

import modelling

for APmodel_name in modelling.model_naming.APmodel_list:

    APsim = modelling.ModelSimController(APmodel_name)
    APsim.simulate()
    myokit.save_state(os.path.join(modelling.RESULT_DIR, 'steady_states',
                                   f'{APmodel_name}_steadystate_APDprot.csv'),
                      APsim.sim.state())

    APsim = modelling.ModelSimController(APmodel_name)
    APsim.sim.set_protocol(myokit.pacing.blocktrain(period=2000,
                                                    duration=0.5,
                                                    offset=50))
    APsim.simulate()
    myokit.save_state(os.path.join(modelling.RESULT_DIR, 'steady_states',
                                   f'{APmodel_name}_steadystate_qNetprot.csv'),
                      APsim.sim.state())

    print(f'{APmodel_name} done!')

for IKr_model in ['Li', 'Lei']:
    sim = modelling.ModelSimController(IKr_model)
    sim.simulate()
    myokit.save_state(os.path.join(modelling.RESULT_DIR, 'steady_states',
                                   f'{IKr_model}_steadystate_Milnes.csv'),
                      sim.sim.state())

    print(f'{IKr_model} done!')
