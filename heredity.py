import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}

def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):
                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)
    
    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")

def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data

def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]

def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """

    list_of_probabilities = list()

    for person in people:

        mother = people[person]["mother"]
        father = people[person]["father"]
        
        if mother in one_gene:
            mother_pass_probability = 0.5
        elif mother in two_genes:
            mother_pass_probability = 1.0 - PROBS["mutation"]
        else:
            mother_pass_probability = PROBS["mutation"]

        father_pass_probability = 0.0
        if father in one_gene:
            father_pass_probability = 0.5
        elif father in two_genes:
            father_pass_probability = 1.0 - PROBS["mutation"]
        else:
            father_pass_probability = PROBS["mutation"]

        if person in two_genes:
            if person in have_trait:
                list_of_probabilities.append(PROBS["trait"][2][True])
            else:
                list_of_probabilities.append(PROBS["trait"][2][False])

            if mother and father:
                list_of_probabilities.append(mother_pass_probability * father_pass_probability)
            else:
                list_of_probabilities.append(PROBS["gene"][2])

        elif person in one_gene:
            if person in have_trait:
                list_of_probabilities.append(PROBS["trait"][1][True])
            else:
                list_of_probabilities.append(PROBS["trait"][1][False])

            if mother and father:
                list_of_probabilities.append((1.0 - mother_pass_probability) * father_pass_probability + mother_pass_probability * (1.0 - father_pass_probability))
            else:
                list_of_probabilities.append(PROBS["gene"][1])
                
        else:
            if person in have_trait:
                list_of_probabilities.append(PROBS["trait"][0][True])
            else:
                list_of_probabilities.append(PROBS["trait"][0][False])

            if mother and father:
                list_of_probabilities.append((1.0 - mother_pass_probability) * (1.0 - father_pass_probability))
            else:
                list_of_probabilities.append(PROBS["gene"][0])

    result = 1.0
    for probability in list_of_probabilities:
        result *= probability

    return result


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """

    for person in probabilities:
        if person in two_genes:
            probabilities[person]["gene"][2] += p
        elif person in one_gene:
            probabilities[person]["gene"][1] += p
        else:
            probabilities[person]["gene"][0] += p

        if person in have_trait:
            probabilities[person]["trait"][True] += p
        else:
            probabilities[person]["trait"][False] += p

def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        
        sum_genes = probabilities[person]["gene"][0] + probabilities[person]["gene"][1] + probabilities[person]["gene"][2]
        if sum_genes > 0.001:
            probabilities[person]["gene"][0] /= sum_genes
            probabilities[person]["gene"][1] /= sum_genes
            probabilities[person]["gene"][2] /= sum_genes
        
        sum_traits = probabilities[person]["trait"][False] + probabilities[person]["trait"][True]
        if sum_traits > 0.001:
            probabilities[person]["trait"][False] /= sum_traits
            probabilities[person]["trait"][True] /= sum_traits

if __name__ == "__main__":
    main()
