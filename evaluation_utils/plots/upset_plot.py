import pandas as pd
from upsetplot import UpSet, from_memberships
import matplotlib.pyplot as plt


def generate_upset_plot(intensities, colname, title, splited=False, 
                        categories=['Center1', 'Center2', 'Center3'],
                        save_plot=False, path_svg=None):
    # Extracting unique features from each center
    features_sets = {}
    for center, details in intensities.items():
        if splited:
            features_sets[center] = set(';'.join(details[colname]).split(';'))
        else:
            features_sets[center] = set(details[colname])

        
    # Preparing data for DataFrame construction
    data = {'value': []}
    for center in categories:
        data[center] = []


    # Combining all unique features from the centers
    all_features = set.union(*features_sets.values())

    # Filling the data dictionary
    for item in all_features:
        data['value'].append(item)
        for center in categories:
            data[center].append(item in features_sets.get(center, []))
    
    # Creating a DataFrame from the data
    df = pd.DataFrame(data)

    # Generating membership list for UpSet plot
    membership_list = df.drop('value', axis=1).astype(bool).apply(lambda row: row.index[row].tolist(), axis=1)
    example = from_memberships(membership_list, data=df['value'])

    # Creating and displaying the UpSet plot
    upset = UpSet(example, subset_size='count', show_counts=True, sort_by='cardinality')
    upset.plot()
    plt.title(title)

    if save_plot:
        plt.savefig(path_svg, dpi=300, format='svg')
        plt.show()
    else:
        plt.show()