import matplotlib.pyplot as plt
import numpy as np
from .plotting.distributions import plot_distance_posteriors, plot_temporal_distribution
from .plotting.skymap import plot_skymap

def plot_association_summary(results_dict, output_file="association_summary.png"):
    '''Create summary plot showing all association statistics'''
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Overlap integrals bar chart
    ax1 = axes[0, 0]
    integrals = ['Spatial\\n(I_Ω)', 'Distance\\n(I_DL)', 'Temporal\\n(I_t)']
    values = [results_dict.get('I_omega', 0),
             results_dict.get('I_dl', 0),
             results_dict.get('I_t', 0)]
    
    bars = ax1.bar(integrals, values, color=['blue', 'green', 'red'])
    ax1.set_ylabel('Overlap Integral Value', fontsize=12)
    ax1.set_title('Individual Overlap Integrals', fontsize=14)
    ax1.set_ylim(0, max(values) * 1.2 if max(values) > 0 else 1)
    
    # Add values on bars
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.3e}', ha='center', va='bottom', fontsize=10)
    
    # Bayes factor visualization
    ax2 = axes[0, 1]
    bf = results_dict.get('bayes_factor', 1.0)
    log_bf = np.log10(bf) if bf > 0 else -10
    
    # Jeffrey's scale for Bayes factors
    ax2.barh([0], [log_bf], color='purple', height=0.5)
    ax2.set_xlabel('log₁₀(Bayes Factor)', fontsize=12)
    ax2.set_title('Bayes Factor', fontsize=14)
    ax2.set_xlim(-5, 5)
    ax2.axvline(0, color='black', linestyle='--', alpha=0.5)
    ax2.set_yticks([])
    
    # Posterior odds
    ax3 = axes[1, 0]
    odds = results_dict.get('posterior_odds', 1.0)
    prob = results_dict.get('confidence', 0.5)
    
    # Pie chart for probability
    sizes = [prob, 1-prob]
    labels = [f'Associated\\n({prob:.1%})', f'Not Associated\\n({(1-prob):.1%})']
    colors = ['green' if prob > 0.5 else 'red', 'lightgray']
    
    ax3.pie(sizes, labels=labels, colors=colors,
           autopct='', startangle=90)
    ax3.set_title(f'Association Probability\\nOdds = {odds:.2e}', fontsize=14)
    
    # Summary statistics table
    ax4 = axes[1, 1]
    ax4.axis('tight')
    ax4.axis('off')
    
    summary_data = [
        ['Metric', 'Value'],
        ['─' * 20, '─' * 20],
        ['Spatial Overlap', f"{results_dict.get('I_omega', 0):.3e}"],
        ['Distance Overlap', f"{results_dict.get('I_dl', 0):.3e}"],
        ['Temporal Overlap', f"{results_dict.get('I_t', 0):.3e}"],
        ['─' * 20, '─' * 20],
        ['Bayes Factor', f"{bf:.3e}"],
        ['Posterior Odds', f"{odds:.3e}"],
        ['Log₁₀ Odds', f"{results_dict.get('log_posterior_odds', -np.inf):.2f}"],
        ['─' * 20, '─' * 20],
        ['P(Associated)', f"{prob:.2%}"],
    ]
    
    table = ax4.table(cellText=summary_data, loc='center',
                     cellLoc='left', colWidths=[0.6, 0.4])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    
    ax4.set_title('Summary Statistics', fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fig

def plot_candidate_ranking(candidates_results, output_file="candidate_ranking.png"):
    '''Plot ranking of multiple EM candidates'''
    
    if not candidates_results:
        return None
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Extract data from results
    n_candidates = len(candidates_results)
    indices = np.arange(n_candidates)
    
    # Sort by posterior odds
    sorted_results = sorted(candidates_results, 
                          key=lambda x: x.get('odds', 0),
                          reverse=True)
    
    # Bar plot of log odds
    ax1 = axes[0]
    log_odds = [np.log10(r['odds']) if r['odds'] > 0 else -10 
                for r in sorted_results]
    colors = ['green' if lo > 0 else 'red' for lo in log_odds]
    
    bars = ax1.barh(indices, log_odds, color=colors, alpha=0.7)
    ax1.set_ylabel('Candidate Rank', fontsize=12)
    ax1.set_xlabel('Log₁₀(Posterior Odds)', fontsize=12)
    ax1.set_title('Candidate Ranking by Posterior Odds', fontsize=14)
    ax1.axvline(0, color='black', linestyle='--', alpha=0.5)
    ax1.set_yticks(indices)
    ax1.set_yticklabels([f'C{i+1}' for i in range(n_candidates)])
    ax1.invert_yaxis()
    
    # Probability plot
    ax2 = axes[1]
    probs = [r.get('probability', 0) for r in sorted_results]
    
    bars = ax2.barh(indices, probs, color=['green' if p > 0.5 else 'red' 
                                           for p in probs], alpha=0.7)
    ax2.set_xlabel('Association Probability', fontsize=12)
    ax2.set_ylabel('Candidate', fontsize=12)
    ax2.set_title('Association Probabilities', fontsize=14)
    ax2.set_yticks(indices)
    ax2.set_yticklabels([f'C{i+1}' for i in range(n_candidates)])
    ax2.invert_yaxis()
    ax2.axvline(0.5, color='black', linestyle='--', alpha=0.5, label='50% threshold')
    ax2.legend()
    ax2.set_xlim(0, 1)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return fig