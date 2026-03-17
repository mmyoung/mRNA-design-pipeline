#!/usr/bin/env python3
"""
Rank Candidates and Generate HTML Report - Enhanced Version
Includes all mRNA design metrics
"""

import sys
import json
from pathlib import Path
from datetime import datetime

# Default weights for scoring - Updated based on Moderna/Pfizer公开参数
# 参考: 
# - N Engl J Med 2020 (Moderna Phase 1)
# - Nature Biotechnology 2021 (mRNA优化策略)
# - Nat Rev Drug Discov 2020 (mRNA疫苗设计)
WEIGHTS = {
    # ========================================
    # 免疫原性 (35%) - 最关键！LNP疫苗最大挑战
    # Moderna: 用 ψ 替代 U, 降低 TLR7/8/9 激活
    # ========================================
    'CpG_frequency':        0.15,   # CpG → TLR9 (最危险)
    'Uridine_content':      0.15,   # Uridine → TLR7/8
    'Immunogenicity_score': 0.05,   # 综合评分
    
    # ========================================
    # 密码子优化 (30%) - 直接影响翻译效率
    # Moderna: CAI > 0.85, 人类最优密码子
    # ========================================
    'CAI':                  0.12,   # 最重要密码子指标
    'tAI':                  0.08,   # tRNA 适配
    'ENC':                  0.05,   # 有效密码子数
    'Codon_pair_bias':      0.05,   # 密码子对
    
    # ========================================
    # GC含量 (15%) - 平衡稳定性和表达
    # Moderna: 50-60% GC
    # ========================================
    'GC_content':           0.10,
    'GC3_content':          0.05,
    
    # ========================================
    # 翻译效率 (10%) - 综合指标
    # ========================================
    'Translation_efficiency': 0.10,
    
    # ========================================
    # 结构和启动 (10%)
    # Moderna: 5'UTR 短, 避免二级结构
    # ========================================
    'MFE_estimate':         0.03,   # 二级结构稳定性
    '5prime_structure_score': 0.04,  # 5'端可及性
    'Kozak_score':          0.03,   # 翻译起始
}

def normalize(value, min_val, max_val, reverse=False):
    """Normalize value to 0-100 range"""
    if max_val == min_val:
        return 50
    norm = (value - min_val) / (max_val - min_val) * 100
    return 100 - norm if reverse else max(0, min(100, norm))

def calculate_composite(cand, ranges):
    """Calculate weighted composite score"""
    score = 0
    total_weight = 0
    
    # Metrics where LOWER is better (reverse)
    reverse_metrics = {'ENC', 'CpG_frequency', 'Uridine_content', 
                      'Immunogenicity_score', 'MFE_estimate'}
    
    for metric, weight in WEIGHTS.items():
        if metric in cand:
            val = cand[metric]
            min_v, max_v = ranges.get(metric, (0, 100))
            rev = metric in reverse_metrics
            norm = normalize(val, min_v, max_v, rev)
            score += norm * weight
            total_weight += weight
    
    return round(score / total_weight, 2) if total_weight > 0 else 0

def get_immunogenicity_badge(score):
    """Get badge color for immunogenicity"""
    if score < 20:
        return 'badge-good', 'Low'
    elif score < 40:
        return 'badge-moderate', 'Moderate'
    else:
        return 'badge-poor', 'High'

def get_rank_class(rank):
    """Get CSS class for ranking"""
    if rank == 1:
        return 'rank-1'
    elif rank == 2:
        return 'rank-2'
    elif rank == 3:
        return 'rank-3'
    return ''

def generate_report(data):
    """Generate HTML report"""
    # Parse input
    if len(sys.argv) < 2:
        print("Usage: rank_and_report.py <scores.json>", file=sys.stderr)
        sys.exit(1)
    
    with open(sys.argv[1], 'r') as f:
        candidates = json.load(f)
    
    if not candidates:
        print("Error: No candidates found", file=sys.stderr)
        sys.exit(1)
    
    # Calculate ranges for normalization
    ranges = {
        'CAI': (min(c['CAI'] for c in candidates), max(c['CAI'] for c in candidates)),
        'tAI': (min(c['tAI'] for c in candidates), max(c['tAI'] for c in candidates)),
        'ENC': (20, 61),
        'GC_content': (30, 70),
        'GC3_content': (30, 80),
        'CpG_frequency': (0, 10),
        'Uridine_content': (15, 40),
        'Immunogenicity_score': (0, 100),
        'MFE_estimate': (min(c['MFE_estimate'] for c in candidates), 
                        max(c['MFE_estimate'] for c in candidates)),
        '5prime_structure_score': (0, 1),
        'Coding_stability': (0, 1),
        'Translation_efficiency': (min(c['Translation_efficiency'] for c in candidates),
                                  max(c['Translation_efficiency'] for c in candidates)),
        'Codon_pair_bias': (min(c['Codon_pair_bias'] for c in candidates),
                            max(c['Codon_pair_bias'] for c in candidates)),
        'Kozak_score': (0, 1),
    }
    
    # Calculate composite scores
    for c in candidates:
        c['composite_score'] = calculate_composite(c, ranges)
    
    # Sort by score (descending)
    candidates.sort(key=lambda x: x['composite_score'], reverse=True)
    
    # Assign ranks
    for i, c in enumerate(candidates):
        c['rank'] = i + 1
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>mRNA Design Pipeline - Enhanced Results</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; 
               background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
               padding: 20px; min-height: 100vh; }}
        .container {{ max-width: 1600px; margin: 0 auto; background: white; 
                     border-radius: 20px; overflow: hidden; box-shadow: 0 20px 60px rgba(0,0,0,0.3); }}
        .header {{ background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%); 
                  color: white; padding: 30px; text-align: center; }}
        .header h1 {{ margin: 0; font-size: 2.5em; }}
        .header p {{ opacity: 0.8; margin-top: 10px; }}
        .stats {{ display: flex; gap: 15px; padding: 20px; background: #f8f9fa; flex-wrap: wrap; justify-content: center; }}
        .stat {{ flex: 1; min-width: 140px; max-width: 200px; text-align: center; padding: 15px; 
                background: white; border-radius: 10px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        .stat .value {{ font-size: 1.8em; font-weight: bold; color: #667eea; }}
        .stat .label {{ color: #666; margin-top: 5px; font-size: 0.9em; }}
        .section {{ padding: 20px; }}
        .section h2 {{ color: #1a1a2e; margin-bottom: 15px; font-size: 1.5em; }}
        .charts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }}
        .chart-box {{ background: #f8f9fa; border-radius: 10px; padding: 15px; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
        th {{ background: #1a1a2e; color: white; padding: 12px 10px; text-align: left; 
             position: sticky; top: 0; white-space: nowrap; }}
        td {{ padding: 10px; border-bottom: 1px solid #eee; }}
        tr:hover {{ background: #f5f5f5; }}
        .rank-1 {{ background: #ffd700 !important; }}
        .rank-2 {{ background: #c0c0c0 !important; }}
        .rank-3 {{ background: #cd7f32 !important; }}
        .badge {{ padding: 3px 8px; border-radius: 10px; font-size: 0.8em; }}
        .badge-good {{ background: #28a745; color: white; }}
        .badge-moderate {{ background: #ffc107; color: #333; }}
        .badge-poor {{ background: #dc3545; color: white; }}
        .score {{ font-weight: bold; color: #667eea; }}
        .metric-section {{ background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 20px; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; }}
        .metric-item {{ background: white; padding: 10px; border-radius: 5px; }}
        .metric-name {{ font-size: 0.85em; color: #666; }}
        .metric-value {{ font-weight: bold; color: #333; }}
        .scroll-table {{ max-height: 600px; overflow-y: auto; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 mRNA Design Pipeline</h1>
            <p>Enhanced Candidate Ranking Report | {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
            <p>Evaluating {len(candidates)} candidate sequences</p>
        </div>
        
        <div class="stats">
            <div class="stat">
                <div class="value">{len(candidates)}</div>
                <div class="label">Total Candidates</div>
            </div>
            <div class="stat">
                <div class="value">{candidates[0]['composite_score']:.1f}</div>
                <div class="label">Top Score</div>
            </div>
            <div class="stat">
                <div class="value">{sum(c['composite_score'] for c in candidates)/len(candidates):.1f}</div>
                <div class="label">Average Score</div>
            </div>
            <div class="stat">
                <div class="value">{max(c['CAI'] for c in candidates):.3f}</div>
                <div class="label">Best CAI</div>
            </div>
            <div class="stat">
                <div class="value">{min(c['Immunogenicity_score'] for c in candidates):.1f}</div>
                <div class="label">Lowest Immunogenicity</div>
            </div>
            <div class="stat">
                <div class="value">{max(c['Translation_efficiency'] for c in candidates):.3f}</div>
                <div class="label">Best Translation</div>
            </div>
        </div>
        
        <div class="section">
            <h2>📊 Score Distribution</h2>
            <div class="charts">
                <div class="chart-box">
                    <div id="barChart"></div>
                </div>
                <div class="chart-box">
                    <div id="histChart"></div>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>📈 Key Metrics Comparison (Top 10)</h2>
            <div class="charts">
                <div class="chart-box">
                    <div id="radarChart"></div>
                </div>
                <div class="chart-box">
                    <div id="caiChart"></div>
                </div>
            </div>
        </div>
        
        <div class="section">
            <h2>🏆 Top 50 Candidates</h2>
            <div class="scroll-table">
                <table>
                    <thead>
                        <tr>
                            <th>Rank</th>
                            <th>Score</th>
                            <th>CAI</th>
                            <th>tAI</th>
                            <th>ENC</th>
                            <th>GC%</th>
                            <th>GC3%</th>
                            <th>CpG%</th>
                            <th>U%</th>
                            <th>Immuno</th>
                            <th>Kozak</th>
                            <th>Trans.Eff</th>
                            <th>MFE</th>
                            <th>5' Str</th>
                            <th>CPS</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    # Add table rows
    for c in candidates[:50]:
        rank_class = get_rank_class(c['rank'])
        immuno_badge, immuno_label = get_immunogenicity_badge(c.get('Immunogenicity_score', 0))
        
        html += f"""                        <tr class="{rank_class}">
                            <td><strong>#{c['rank']}</strong></td>
                            <td class="score">{c['composite_score']:.2f}</td>
                            <td>{c['CAI']:.4f}</td>
                            <td>{c['tAI']:.4f}</td>
                            <td>{c['ENC']:.1f}</td>
                            <td>{c['GC_content']:.1f}%</td>
                            <td>{c['GC3_content']:.1f}%</td>
                            <td>{c['CpG_frequency']:.2f}%</td>
                            <td>{c['Uridine_content']:.1f}%</td>
                            <td><span class="badge {immuno_badge}">{c.get('Immunogenicity_score', 0):.1f}</span></td>
                            <td>{c['Kozak_score']:.1f}</td>
                            <td>{c['Translation_efficiency']:.3f}</td>
                            <td>{c['MFE_estimate']:.1f}</td>
                            <td>{c['5prime_structure_score']:.2f}</td>
                            <td>{c['Codon_pair_bias']:.1f}</td>
                        </tr>
"""
    
    # Get data for charts
    top10 = candidates[:10]
    
    html += f"""                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="section metric-section">
            <h2>📋 Scoring Weights</h2>
            <div class="metric-grid">
                <div class="metric-item"><div class="metric-name">CAI (Codon Adaptation)</div><div class="metric-value">10%</div></div>
                <div class="metric-item"><div class="metric-name">tAI (tRNA Adaptation)</div><div class="metric-value">10%</div></div>
                <div class="metric-item"><div class="metric-name">ENC</div><div class="metric-value">5%</div></div>
                <div class="metric-item"><div class="metric-name">Codon Pair Bias</div><div class="metric-value">5%</div></div>
                <div class="metric-item"><div class="metric-name">GC Content</div><div class="metric-value">10%</div></div>
                <div class="metric-item"><div class="metric-name">GC3</div><div class="metric-value">5%</div></div>
                <div class="metric-item"><div class="metric-name">CpG Frequency</div><div class="metric-value">10%</div></div>
                <div class="metric-item"><div class="metric-name">Uridine Content</div><div class="metric-value">10%</div></div>
                <div class="metric-item"><div class="metric-name">Immunogenicity</div><div class="metric-value">5%</div></div>
                <div class="metric-item"><div class="metric-name">Translation Efficiency</div><div class="metric-value">10%</div></div>
                <div class="metric-item"><div class="metric-name">5' Structure</div><div class="metric-value">5%</div></div>
                <div class="metric-item"><div class="metric-name">Kozak</div><div class="metric-value">5%</div></div>
            </div>
        </div>
    </div>
    
    <script>
        // Bar chart - Top 10 scores
        new Plotly.newPlot('barChart', [{{
            x: {json.dumps([f"#{c['rank']}" for c in top10])},
            y: {json.dumps([c['composite_score'] for c in top10])},
            type: 'bar',
            marker: {{ color: '#667eea' }}
        }}], {{ 
            title: 'Top 10 - Composite Scores',
            yaxis: {{ title: 'Score' }},
            margin: {{ t: 40, b: 40, l: 50, r: 20 }}
        }});
        
        // Histogram - Score distribution
        new Plotly.newPlot('histChart', [{{
            x: {json.dumps([c['composite_score'] for c in candidates])},
            type: 'histogram',
            marker: {{ color: '#764ba2' }}
        }}], {{
            title: 'Score Distribution',
            xaxis: {{ title: 'Composite Score' }},
            yaxis: {{ title: 'Count' }},
            margin: {{ t: 40, b: 40, l: 50, r: 20 }}
        }});
        
        // CAI comparison
        new Plotly.newPlot('caiChart', [{{
            x: {json.dumps([f"#{c['rank']}" for c in top10])},
            y: {json.dumps([c['CAI'] for c in top10])},
            type: 'bar',
            marker: {{ color: '#28a745' }}
        }}], {{
            title: 'Top 10 - CAI Values',
            yaxis: {{ title: 'CAI', range: [0, 0.35] }},
            margin: {{ t: 40, b: 40, l: 50, r: 20 }}
        }});
        
        // Radar chart - normalized metrics for top 3
        var top3 = {json.dumps([c['CAI'] for c in top10[:3]])};
        var metrics = ['CAI', 'tAI', 'GC%', 'CpG%', 'U%', 'Kozak', 'TransEff'];
        
        new Plotly.newPlot('radarChart', [
"""
    
    # Add traces for top 3
    colors = ['#ffd700', '#c0c0c0', '#cd7f32']
    for i, c in enumerate(top10[:3]):
        html += f"""            {{
            r: [{c['CAI']*300}, {c['tAI']*300}, {c['GC_content']/100*100}, {c['CpG_frequency']*10}, 
                {c['Uridine_content']*2.5}, {c['Kozak_score']*100}, {c['Translation_efficiency']*100}],
            theta: metrics,
            fill: 'toself',
            name: '#{c["rank"]}',
            marker: {{ color: '{colors[i]}' }}
        }},
"""
    
    html += """        ], {
            title: 'Top 3 - Metrics Comparison',
            polar: { radialaxis: { visible: true, range: [0, 100] } },
            margin: { t: 40, b: 40, l: 50, r: 20 }
        });
    </script>
</body>
</html>"""
    
    print(html)

if __name__ == '__main__':
    generate_report(None)
