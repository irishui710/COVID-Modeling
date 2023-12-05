population = 325_000_000  # 325 million people in the US

# Define the δ values to study
delta_values = [0.24, 0.22, 0.20, 0.18, 0.16]

# Initialize lists to store the results for each δ value
max_hospitalization_list = []
day_max_hospitalization_list = []
end_outbreak_day_list = []
percentage_infected_list = []
total_deaths_list = []

# Set pc to 0.3 for this scenario
pc_ = 0.3

# Below: My Code
# Define the values of δ to study
transmission_rates = [0.24, 0.22, 0.2, 0.18, 0.16]

# Create lists to store results for each δ
cumulative_deaths_list = []
daily_hospitalizations_list = []
susceptible_percentage_list = []

# Iterate through different δ values
for delta_ in transmission_rates:
    # Run the SEAIHDR model for the current δ
    S, E, A, I, H, R, D, _ = covid19_model(delta_, pc_, pca_, R0_, h0_, r0_)

    # Calculate daily hospitalizations per 100K people
    population = 325_000_000  # 325 million people in the US
    daily_hospitalizations = H * 100_000  # Scale to per 100K people

    # Calculate cumulative deaths
    cumulative_deaths = D * population

    # Calculate the percentage of the population still susceptible to COVID-19
    susceptible_percentage = S

    # Store the results for the current δ
    cumulative_deaths_list.append(cumulative_deaths)
    daily_hospitalizations_list.append(daily_hospitalizations)
    susceptible_percentage_list.append(susceptible_percentage)

# Create plots for each parameter
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot cumulative deaths
axes[0].set_title('Cumulative Deaths')
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Cumulative Deaths')
for i, delta in enumerate(transmission_rates):
    axes[0].plot(cumulative_deaths_list[i], label=f'delta = {delta}')
axes[0].legend()

# Plot daily hospitalizations per 100K people
axes[1].set_title('Daily Hospitalizations')
axes[1].set_xlabel('Days')
axes[1].set_ylabel('Daily Hospitalizations per 100K People')
for i, delta in enumerate(transmission_rates):
    axes[1].plot(daily_hospitalizations_list[i], label=f'delta = {delta}')
axes[1].legend()

# Plot the percentage of the population still susceptible
axes[2].set_title('Susceptible Population')
axes[2].set_xlabel('Days')
axes[2].set_ylabel('Percentage of Population Susceptible')
for i, delta in enumerate(transmission_rates):
    axes[2].plot(susceptible_percentage_list[i], label=f'delta = {delta}')
axes[2].legend()

plt.tight_layout()
plt.show()

# Iterate through different δ values
for delta_ in delta_values:
    # Run the SEAIHDR model for the current δ
    S, E, A, I, H, R, D, _ = covid19_model(delta_, pc_, pca_, R0_, h0_, r0_)

    # Calculate daily new cases per 100K people
    population = 325_000_000  # 325 million people in the US
    daily_new_cases = -np.diff(S) * population
    daily_new_cases /= 100_000  # Scale to per 100K people

    # Calculate the outcomes
    max_hospitalization = np.max(H) * 100_000
    day_max_hospitalization = np.argmax(H)
    total_infected = A + I + H
    end_day = np.min(np.where(total_infected[1:] * 100_000 < 100))
    percentage_infected = 1 - S[-1]
    total_deaths = D[-1] * population

    # Append the results to the respective lists
    max_hospitalization_list.append(max_hospitalization)
    day_max_hospitalization_list.append(day_max_hospitalization)
    end_outbreak_day_list.append(end_day)
    percentage_infected_list.append(percentage_infected)
    total_deaths_list.append(total_deaths)

# Print the results for each δ value
for i, delta in enumerate(delta_values):
    print(f"Results for delta = {delta}:")
    print(f"i. Maximum Hospitalization per 100K People: {max_hospitalization_list[i]:.2f}")
    print(f"ii. Day of Maximum Hospitalization: {day_max_hospitalization_list[i]}")
    print(f"iii. Day of End of Outbreak: {end_outbreak_day_list[i]}")
    print(f"iv. Percentage of Population Infected: {percentage_infected_list[i]:.2f}")
    print(f"v. Total Deaths: {total_deaths_list[i]:.0f}\n")


# PART D - Trial and Error
import numpy as np

# Function to run the SEAIHDR model and return total deaths
def run_simulation(delta):
    # Run the SEAIHDR simulation with the given delta
    S, E, A, I, H, R, D, _ = covid19_model(delta, pc_, pca_, R0_, h0_, r0_)
    return D[-1] * population

# Define the target range for total deaths
target_lower = 100000
target_upper = 200000

# Initialize variables for delta range
delta_lower = 0.0
delta_upper = 1.0

# Specify the precision for the estimate
precision = 0.001

# Initialize variables
best_delta = 0
best_deaths = float('inf')  # Initialize with a large value

# Find the lower bound of delta
while delta_lower < 1.0:
    deaths_lower = run_simulation(delta_lower)
    if deaths_lower < target_lower:
        delta_lower += precision
    else:
        break

# Find the upper bound of delta
while delta_upper > 0.0:
    deaths_upper = run_simulation(delta_upper)
    if deaths_upper > target_upper:
        delta_upper -= precision
    else:
        break

# Print the delta range
print(f"Lower Bound of Delta: {delta_lower:.3f}")
print(f"Upper Bound of Delta: {delta_upper:.3f}")



