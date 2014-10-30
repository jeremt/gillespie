import random
import math


def algorithm(molecules, system, constants, max_time=100., normalized=True):
    """
    Simulates chemical reactions using Gillespie algorithm.
    :param molecules: The number of molecules of each type (as a dictionary)
    :param system: A list of dictionaries which contain the propensity equation and the reaction to apply
    :param constants: All the constants involved in the algorithm
    """

    # Init default values
    molecules_0 = molecules.copy()
    variables = constants.copy()
    times = []
    points = dict((name, []) for name in molecules)
    current_time = 0
    while sum(molecules.values()) > 1 and current_time < max_time:

        # Save the points and the current time
        for name in molecules:
            points[name].append(molecules[name] / (sum(molecules.values()) if normalized else 1))
        times.append(current_time)

        # Compute the propensities of each reaction
        variables.update(molecules)
        propensities = [eval(eq['propensity'], variables.copy()) for eq in system]
        if sum(propensities) == 0:
            break

        # Randomly choose the reaction to apply
        reaction = random.uniform(0.0, sum(propensities))
        for i, eq in enumerate(system):
            p = sum(propensities[0:i+1])
            if p > 0 and reaction < p:
                for name in eq['reaction']:
                    molecules[name] += eq['reaction'][name]

        # Increase the time
        current_time += -(1. / sum(propensities)) * math.log(random.random())

    return {
        'molecules': molecules_0,
        'points': points,
        'times': times,
        'constants': constants,
    }
