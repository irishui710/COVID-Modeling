import numpy as np
import matplotlib.pyplot as plt

# Suggested default values
global delta_, pc_, pca_, R0_, h0_, r0_
delta_ = 1
pc_ = 0.1
pca_ = 0
R0_ = 5.7
h0_ = 0.1  # Toggled between 1 and 0.1 for problem 1, part a
r0_ = 0

# Common Data
global pa_, ph_, mh_, fa_, fc_, fh_, te_, ta_, ti_, tc_, th_, target_, maxdays_
pa_ = 0.4
ph_ = 0.016
mh_ = 0.25
fa_ = 0.75
fc_ = 0.1
fh_ = 0.05
te_ = 5
ta_ = 8
ti_ = 10
tc_ = 3
th_ = 10
target_ = 60
maxdays_ = 1000

## Derived params
global eta_, alpha_, gamma_, nu_, sigma_, beta_
eta_ = 1 / te_
alpha_ = 1 / ta_
gamma_ = 1 / ti_
nu_ = 1 / th_
sigma_ = 1 / tc_
beta_ = R0_ / (fa_ * pa_ * ta_ + (1 - pa_) * ta_)


## FUNCTION FOR THE DIFFERENTIAL EQUATION
def yprime(y):
    global eta_, alpha_, gamma_, nu_, sigma_, beta_

    # split out components
    S = y[0]
    E = y[1]
    A = y[2]
    I1 = y[3]
    I2 = y[4]
    I = I1 + I2
    H = y[5]
    X = delta_ * (fa_ * (1 - pca_) * A + (1 - pc_) * I) + fc_ * (pca_ * A + pc_ * I) + fh_ * H

    # compute derivatives
    Sp = -beta_ * X * S
    Ep = -Sp - eta_ * E
    Ap = pa_ * eta_ * E - alpha_ * A
    I1p = (1 - pa_ - ph_) * eta_ * E - gamma_ * I1
    I2p = ph_ * eta_ * E - sigma_ * I2
    Hp = sigma_ * I2 - nu_ * H
    Rp = alpha_ * A + gamma_ * I1 + (1 - mh_) * nu_ * H
    Dp = mh_ * nu_ * H

    # assemble derivative
    yp = np.array([Sp, Ep, Ap, I1p, I2p, Hp, Rp, Dp])

    return yp


## FUNCTION FOR rk4
def rk4(dt, y0):
    # y0 is a column vector of initial conditions at t
    # y is a column vector of values at t+dt
    k1 = yprime(y0)
    k2 = yprime(y0 + 0.5 * dt * k1)
    k3 = yprime(y0 + 0.5 * dt * k2)
    k4 = yprime(y0 + dt * k3)
    y = y0 + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return y


def covid19_model(delta_in, pc_in, pca_in, R0_in, h0_in, r0_in):
    #
    # returns [S,E,A,I,H,R,D,lambdaa] = covid19_sim(delta_in,pc_in,pca_in,R0_in,h0_in,r0_in)
    #
    # runs a simulation of an SEAIHDR model
    #
    # S: susceptible
    # E: exposed
    # A: asymptomatic
    # I: infective (symptomatic)
    # H: hospitalized
    # R: recovered
    # D: deceased
    # lambdaa: initial exponential growth rate
    #
    # delta is the fractional contact rate for social distancing
    # pa is the fraction of infectives who are asymptomatic
    # pc is the fraction of symptomatics who are tested
    # pca is the fraction of asymptomatics who are tested
    # R0 is the basic reproductive number with no detection or intervention
    # ph is the fraction of infectives who require hospitalization
    # mh is the fraction of hospitalized patients who die
    # fa_, fc_, fh_ are the infectivities of A, C, and H, relative to that of I
    # te, ta, ti, tc_, th are the expected durations of classes E, A, I1, I2, H
    # h0 is the number of initial hospitalized infectives per hundred thousand
    # r0 is the fraction of the population that is initially immune
    # target_ is the value of A+I+H per hundred thousand defined as the end condition
    #
    # result is a matrix: columns are SEAI1I2HRD, rows are 0:days
    # Outputs S, E, A, I, H, R, D column vectors (rows are values per day) and scalar lambda
    #
    # adapted from UNL

    ## DATA

    # update data
    global pa_, ph_, mh_, fa_, fc_, fh_, te_, ta_, ti_, tc_, th_, target_, maxdays_
    global delta_, pc_, pca_, R0, h0_, r0_
    delta_ = delta_in
    pc_ = pc_in
    pca_ = pca_in
    R0_ = R0_in
    h0_ = h0_in
    r0_ = r0_in

    ## INITIALIZATION
    global eta_, alpha_, gamma_, nu_, sigma_, beta_

    # set derived parameters
    h0_rate = h0_ / 100000
    target_rate = target_ / 100000

    eta_ = 1 / te_
    alpha_ = 1 / ta_
    gamma_ = 1 / ti_
    nu_ = 1 / th_
    sigma_ = 1 / tc_
    beta_ = R0_ / (fa_ * pa_ * ta_ + (1 - pa_) * ta_)

    M = np.array([[-eta_, fa_ * beta_, beta_, beta_, fh_ * beta_], [pa_ * eta_, -alpha_, 0, 0, 0],
                  [(1 - pa_ - ph_) * eta_, 0, -gamma_, 0, 0], [ph_ * eta_, 0, 0, -sigma_, 0], [0, 0, 0, sigma_, -nu_]])
    eval, evec = np.linalg.eig(M)
    lambdaa = np.amax(eval)

    i20 = h0_rate * (lambdaa + nu_) / sigma_
    e0 = i20 * (lambdaa + sigma_) / (ph_ * eta_)
    a0 = e0 * pa_ * eta_ / (lambdaa + alpha_)
    i10 = e0 * (1 - pa_ - ph_) * eta_ / (lambdaa + alpha_)

    # set up results data structure with Y=[S,E,A,I1,I2,H,R,D]

    results = np.zeros((maxdays_ + 1, 8))
    Y = np.array([0, e0, a0, i10, i20, h0_rate, r0_, 0])
    Y[0] = 1 - np.sum(Y)
    results[0, :] = Y

    y = np.transpose(Y)
    summ = e0 + a0 + i10 + i20 + h0_rate  # summ is used to help determine the end condition

    ## COMPUTATION

    for t in range(maxdays_):
        # y is a column vector, Y^T
        y = rk4(1, y)
        Y = np.transpose(y)
        results[t + 1, :] = Y
        if np.sum(Y[1:5]) > np.minimum(target_rate, summ):
            summ = np.sum(Y[1:5])
        else:
            results = results[0:t, :]
            break

    S = results[:, 0]
    E = results[:, 1]
    A = results[:, 2]
    I = results[:, 3] + results[:, 4]
    H = results[:, 5]
    R = results[:, 6]
    D = results[:, 7]

    return S, E, A, I, H, R, D, lambdaa


'''Plotting'''
# Set the parameter values
pc_ = 0 # set pc to 0

# Run the SEAIHDR model
S, E, A, I, H, R, D, _ = covid19_model(delta_, pc_, pca_, R0_, h0_, r0_)

# Calculate daily new cases per 100K people
population = 325_000_000  # 325 million people in the US
daily_new_cases = -np.diff(S) * 100_000
#daily_new_cases /= 100_000  # Scale to per 100K people

# Calculate daily hospitalizations per 100K people
daily_hospitalizations = H * 100_000  # Scale to per 100K people

# Create the required plots
days = np.arange(len(S))

# Plot fractional values of S, E, A + I + H, R, and D
plt.figure(figsize=(12, 6))
plt.plot(days, S, label='S')
plt.plot(days, E, label='E')
plt.plot(days, A + I + H, label='A + I + H')
plt.plot(days, R, label='R')
plt.plot(days, D, label='D')
plt.xlabel('Days')
plt.ylabel('Fraction of Population')
plt.title('Fractional Values of SEAIHDR Model with pc = 0')
plt.legend()
plt.grid(True)
plt.show()

# Plot cumulative deaths
plt.figure(figsize=(12, 6))
plt.plot(days, D * population)  # Scale to cumulative deaths per 100K people
plt.xlabel('Days')
plt.ylabel('Cumulative Deaths')
plt.title('Cumulative Deaths')
plt.grid(True)
plt.show()

# Plot daily new cases per 100K people
plt.figure(figsize=(12, 6))
plt.plot(days[1:], daily_new_cases)
plt.xlabel('Days')
plt.ylabel('Daily New Cases per 100K People')
plt.title('Daily New Cases')
plt.grid(True)
plt.show()

# Plot daily hospitalizations per 100K people
plt.figure(figsize=(12, 6))
plt.plot(days, daily_hospitalizations)
plt.xlabel('Days')
plt.ylabel('Daily Hospitalizations per 100K People')
plt.title('Daily Hospitalizations')
plt.grid(True)
plt.show()

# Print Model Estimates
# Calculate and print model estimates
ph_max = np.amax(daily_hospitalizations)
day_max_hospitalization = np.argmax(H) + 1
total_infected = A + I + H
end_period = total_infected[-100:-1]
end_day = np.argmax(end_period < (100 / 100_000)) + len(total_infected) - 99
percentage_infected = 1-S[-1]
total_deaths = D[-1] * population

print(f"Maximum Hospitalization Requirement per 100K People: {ph_max:.2f}")
print(f"Day of Maximum Hospitalization Requirement: {day_max_hospitalization}")
print(f"Day When Total Infected Drops Below 100 per 100K: {end_day}")
print(f"Percentage of Population Infected: {percentage_infected:.4f}")
print(f"Total Population Deaths: {total_deaths:.0f}")



