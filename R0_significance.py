# Define the values of R0 to study
R0_values = [3, 4, 5, 6]

# Create lists to store results for each R0
cumulative_deaths_list = []
daily_hospitalizations_list = []
susceptible_percentage_list = []

# Iterate through different R0 values
for R0_ in R0_values:
    # Run the SEAIHDR model for the current R0
    S, E, A, I, H, R, D, _ = covid19_model(delta_, pc_, pca_, R0_, h0_, r0_)

    # Calculate daily hospitalizations per 100K people
    population = 325_000_000  # 325 million people in the US
    daily_hospitalizations = H * 100_000  # Scale to per 100K people

    # Calculate cumulative deaths
    cumulative_deaths = D * population

    # Calculate the percentage of the population still susceptible to COVID-19
    susceptible_percentage = S

    # Store the results for the current R0
    cumulative_deaths_list.append(cumulative_deaths)
    daily_hospitalizations_list.append(daily_hospitalizations)
    susceptible_percentage_list.append(susceptible_percentage)

# Create plots for each parameter
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot cumulative deaths
axes[0].set_title('Cumulative Deaths')
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Cumulative Deaths')
for i, R0 in enumerate(R0_values):
    axes[0].plot(cumulative_deaths_list[i], label=f'R0 = {R0}')
axes[0].legend()

# Plot daily hospitalizations per 100K people
axes[1].set_title('Daily Hospitalizations')
axes[1].set_xlabel('Days')
axes[1].set_ylabel('Daily Hospitalizations per 100K People')
for i, R0 in enumerate(R0_values):
    axes[1].plot(daily_hospitalizations_list[i], label=f'R0 = {R0}')
axes[1].legend()

# Plot the percentage of the population still susceptible
axes[2].set_title('Susceptible Population')
axes[2].set_xlabel('Days')
axes[2].set_ylabel('Percentage of Population Susceptible')
for i, R0 in enumerate(R0_values):
    axes[2].plot(susceptible_percentage_list[i], label=f'R0 = {R0}')
axes[2].legend()

plt.tight_layout()
plt.show()


def doubling_time(hospitalizations):
    # Find the day at which hospitalizations double
    for day, hosp_count in enumerate(hospitalizations):
        if day > 0 and hosp_count >= 2 * hospitalizations[0]:
            return day + 1


# Iterate through different R0 values
for R0_ in R0_values:
    # Run the SEAIHDR model for the current R0
    S, E, A, I, H, R, D, _ = covid19_model(delta_, pc_, pca_, R0_, h0_, r0_)

    # Calculate daily hospitalizations per 100K people
    population = 325_000_000  # 325 million people in the US
    daily_hospitalizations = H * 100_000  # Scale to per 100K people

    # Print the doubling time for daily hospitalizations
    doubling_day = doubling_time(daily_hospitalizations)
    print(f'Doubling time for R0={R0_}: {doubling_day} days')