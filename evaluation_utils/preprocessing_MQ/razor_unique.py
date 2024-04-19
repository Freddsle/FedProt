import pandas as pd
import networkx as nx


def find_leading(features, razor_uniq_dict, more_then_half=False):
    # Filter the counts for the genes in the set
    relevant_counts = {feature: razor_uniq_dict[feature] for feature in features}
    # Find the max count
    max_count = max(relevant_counts.values())
    
    if more_then_half:
        leading_genes = [feature for feature, count in relevant_counts.items() if count >= max_count / 2]
    
    else:
        leading_genes = [feature for feature, count in relevant_counts.items() if count  == max_count]

    return set(sorted(leading_genes))


def build_peptide_gene_graph(peptides_mapping_dict):
    """
    Build a graph where peptides and proteins are nodes.
    """
    G = nx.Graph()
    for peptide, features in peptides_mapping_dict.items():
        for feature in features:
            G.add_node(peptide, type='peptide')
            G.add_node(feature, type='feature')
            G.add_edge(peptide, feature)
    return G


def sort_by_count(features_list, razor_uniq_dict):
    if len(features_list) == 1:
        return features_list
    # sort by razor unique count, from biggest to smallest
    return sorted(list(set(features_list)), key=lambda x: razor_uniq_dict[x], reverse=True)


def peptide_grouping(merged_mapping, feature_column):

    df_exploded = merged_mapping.assign(**{feature_column: merged_mapping[feature_column].str.split(';')}).explode(feature_column)

    unique_razor = df_exploded[feature_column].value_counts().rename_axis(feature_column).reset_index(name='Unique_razor')
    unique_counts = merged_mapping[merged_mapping[feature_column].str.contains(';') == False][feature_column].value_counts().rename_axis(feature_column).reset_index(name='Unique')

    result = pd.merge(unique_razor, unique_counts, on=feature_column, how='left').fillna({'Unique': 0})

    peptides_mappings = pd.Series(merged_mapping[feature_column].values, index=merged_mapping['Sequence']).to_dict()
    peptides_mapping_dict = {key: set(value.split(';')) for key, value in peptides_mappings.items()}

    razor_uniq_dict = pd.Series(result['Unique_razor'].values, index=result[feature_column]).to_dict()
    unique_genes = result[result['Unique'] > 0][feature_column].tolist()

    G = build_peptide_gene_graph(peptides_mapping_dict)
    connected_components = nx.connected_components(G)

    final_protein_groups = []

    for component in connected_components:
        component_copy = set(component)  
        features = [node for node in component_copy if G.nodes[node]['type'] == 'feature']
        peptides = [node for node in component_copy if G.nodes[node]['type'] == 'peptide']
        
        # find proteins with max count
        for feature in features:
            
            if len(peptides) == 0 or len(features) == 0:
                print('Error')
                print(f'Peptides: {peptides}')
                print(f'Features: {features}')
                raise ValueError('Peptides or features are not empty')
        
            leading = list(find_leading(features, razor_uniq_dict))
            leading_unique = set(leading) & set(unique_genes)
            if len(leading) > 1 and  len(leading_unique) > 0:
                # if intersect with unique_genes - take the first unique gene
                leading = list(leading_unique)[0]
            else:
                leading = leading[0]
                    
            
            # Get all peptides directly connected to the leading protein
            leading_peptides = [peptide for peptide in peptides if G.has_edge(leading, peptide)]
            # Find other proteins connected only to these leading peptides
            other = set()

            for peptide in leading_peptides:
                connected_proteins = set(G.neighbors(peptide)) & set(features)  # Proteins connected to this peptide
                
                # Filter out proteins that are connected to peptides not in leading_peptides
                valid_features = set()
                for protein in connected_proteins:
                    protein_peptides = set(G.neighbors(protein)) & set(peptides)  # Peptides connected to this protein
                    if protein_peptides.issubset(set(leading_peptides)):
                        valid_features.add(protein)

                # Update 'other' with proteins that meet the criteria
                other.update(valid_features)

            other.discard(leading)
            leading_group = list(other) + [leading]

            razor_feature = find_leading(leading_group, razor_uniq_dict, more_then_half=False)
            if len(razor_feature) > 1:
                # keep only the first unique if there are more than one and unique is present
                if len(razor_feature & set(unique_genes)) > 0:
                    razor_feature = sort_by_count(list(razor_feature & set(unique_genes)), razor_uniq_dict)
            razor_feature = list(razor_feature)[0]
            
            final_protein_group = {
                'features': sort_by_count(leading_group, razor_uniq_dict),
                'peptides': sorted(list(set(leading_peptides))),
                'razor_feature': razor_feature,
                'major_features': sort_by_count(list(find_leading(leading_group, razor_uniq_dict, more_then_half=True)), razor_uniq_dict),
            }
            
            final_protein_groups.append(final_protein_group)
                
            # Remove the proteins and peptides from the component_copy
            to_remove = set(leading_group) | set(leading_peptides)
            component_copy -= to_remove
            features = [feature for feature in features if feature not in set(leading_group)]
            peptides = [peptide for peptide in peptides if peptide not in set(leading_peptides)]

            if len(peptides) == 0 and len(features) == 0:
                break
        
        if len(peptides) == 0 and len(features) == 0:
            continue

    return final_protein_groups