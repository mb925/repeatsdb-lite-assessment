import pandas as pd
import config as cfg
import matplotlib.pyplot as plt

def plot_evaluation():
    df = pd.read_csv(cfg.data['evaluation'] + '/evaluation.csv', sep=',')
    df = df.loc[(df['unit_accuracy'] != 'No matches') & (df['unit_accuracy'] != 'No uniprot for this PDB')]
    df = df.astype({'unit_accuracy': 'float', 'unit_precision': 'float','unit_recall': 'float',
                    'region_accuracy': 'float','region_precision': 'float','region_recall': 'float',
                    'classe': 'float', 'topology': 'float', 'fold': 'float', 'clan': 'float'})

    # df[['unit_accuracy','unit_precision','unit_recall']].plot(kind='box')
    # plt.savefig(cfg.data['plots'] + "/unit_phase.png")
    # df[['region_accuracy','region_precision','region_recall']].plot(kind='box')
    # plt.savefig(cfg.data['plots'] + "/region_phase.png")
    df[['classe']].value_counts().sort_values().plot(kind='bar')

    plt.savefig(cfg.data['plots'] + "/class_evaluation.png")

    df[['topology']].value_counts().sort_values().plot(kind='bar')
    plt.savefig(cfg.data['plots'] + "/topology_evaluation.png")

    # it is putting the less amount first
    df[['fold']].value_counts().sort_values().plot(kind='bar')
    plt.savefig(cfg.data['plots'] + "/fold_evaluation.png")

    df[['clan']].value_counts().sort_values().plot(kind='bar')
    plt.savefig(cfg.data['plots'] + "/clan_evaluation.png")

if __name__ == '__main__':
    plot_evaluation()