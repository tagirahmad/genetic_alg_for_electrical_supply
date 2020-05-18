from alg import *


def expenses_calculation(к_л, к_эа, к_тп, с_ал, с_пт, с_пл, с_аэа, с_атп):
    expenses = (
        0.15 * (к_л + к_эа + к_тп)
        + с_ал
        + с_пт
        + с_пл
        + с_аэа
        + с_атп
    )
    return expenses


workshop_container = WorkshopsContainer()

result = expenses_calculation(
    к_л=workshop_container.к_л,
    к_эа=workshop_container.к_эа,
    к_тп=workshop_container.к_тп, 
    с_ал=workshop_container.с_ал,
    с_пт=workshop_container.с_пт,
    с_пл=workshop_container.с_пл,
    с_аэа=workshop_container.с_аэа,
    с_атп=workshop_container.с_атп
)

print(result, 'result')


import pygad
import numpy

num_generations = 50 # Number of generations.
sol_per_pop = 8 # Number of solutions in the population.
num_parents_mating = 4 # Number of solutions to be selected as parents in the mating pool.

# Parameters of the mutation operation.
mutation_percent_genes = 10 # Percentage of genes to mutate. This parameter has no action if the parameter mutation_num_genes exists.
mutation_num_genes = None # Number of genes to mutate. If the parameter mutation_num_genes exists, then no need for the parameter mutation_percent_genes.

parent_selection_type = "tournament" # Type of parent selection.

crossover_type = "two_points" # Type of the crossover operator.

mutation_type = "scramble" # Type of the mutation operator.

keep_parents = 1 # Number of parents to keep in the next population. -1 means keep all parents and 0 means keep nothing.

init_range_low = -2
init_range_high = 5


