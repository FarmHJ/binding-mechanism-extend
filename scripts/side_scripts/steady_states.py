import myokit
import os

import modelling

for APmodel_name in modelling.model_naming.APmodel_list:
    # fname = f'{APmodel_name}_steadystate_APDprot.csv'
    fname = f'{APmodel_name}_steadystate_qNetprot.csv'
# for APmodel_name in ['Li', 'Lei']:
#     fname = f'{APmodel_name}_steadystate_Milnes.csv'

    APsim = modelling.ModelSimController(APmodel_name)
    APsim.sim.set_protocol(myokit.pacing.blocktrain(period=2000,
                                                    duration=0.5,
                                                    offset=50))
    APsim.simulate()
    myokit.save_state(os.path.join(modelling.RESULT_DIR, 'steady_states',
                                   fname),
                      APsim.sim.state())
    print(f'{APmodel_name} done!')
