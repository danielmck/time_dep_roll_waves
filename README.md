# Time Dependent Roll Waves
A number of drivers for solving evoling waveforms in debris flow models using a hyperbolic solver. Two distinct cases are modelled, those of waves on a slope with constant angle in a periodic box and a case of inflow into a long domain. The inflow cases are contained in directories ending '_inflow'.

The directories including 'channel_roll_wave' consist of a two equation model considering flow height and velocity.
The four equation model is in directories named 'four_eqn_var_rho' and considers heigh and velocity as well as basal fluid pressure and particle volume fraction. 
Other directories starting with var_rho include variations of the four equation model used for testing different aspects. The '_vis' suffix indicates that an in plane stress "viscous" term has been added in an attempt to dampen down instabilities present.
