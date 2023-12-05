# COVID-Modeling
SEAIHDR Model for COVID pandemic, can be adapted to predict outcomes of any infectious disease outbreak

Modeling infectious diseases, initial conditions: 

The variables of our SEAIHDR model are as follows:
• S: number of susceptible people
• E: number of exposed people
• A: number of asymptomatic people
• I: number of infected (symptomatic) people • H: number of hospitalized people
• R: number of recovered people
• D: number of deceased people

The parameters in the model are as follows:
• δ is the fractional contact/transmission rate (varies from 0 to 1 to account for mitigation strategies)
• pa is the fraction of infected people who are asymptomatic
• pc is the fraction of symptomatic people who are tested
• pca is the fraction of asymptomatic peple who are tested
• R0 is the basic reproductive number with no detection or intervention
• ph is the fraction of infected people who require hospitalization
• mh is the fraction of hospitalized patients who die
• fa,fc,fh are the infectivities of A, C, and H, relative to that of I
• te, ta, ti, tc, th are the expected durations of classes E, A, I1, I2, H
• h0 is the number of initial hospitalized infectives per hundred thousand
• r0 is the fraction of the population that is initially immune
• target is the value of A+I+H per hundred thousand defined as the end condition of the pandemic


covid19_model_all.py includes the SEAIHDR model with numerical method of 4th order Runge Kutta. It also contains code to plot fractional values of SEAIHDR, cumulative deaths, and maximum hospitalizations.

Flattening_the_curve.py includes code to explore the effect of changing transmission rates on all dependent variables explored above. 

blunting_the_surge.py includes code that manipulate parameters in the SEAIHDR to explore the circumstances where we could blunt the surge and effectively prevent the outbreak from becoming an epidemic/pandemic.

R0_significance explores the effect of the R0 value on the outcome of the pandemic and underscores the importance of calculating and monitoring R0 during a pandemic.
