#!/usr/bin/env python3
"""
Step 5: Generate HTML Report
Create interactive HTML report with visualizations
"""

import sys
import json
import base64
from pathlib import Path
from datetime import datetime

def create_plotly_figure(data, chart_type='bar'):
    """Create Plotly figure as base64 image"""
    # This is a placeholder - actual implementation would use plotly
    # For now, we'll create a simple HTML-based chart
    return ""

def generate_radar_chart(candidate, metrics):
    """Generate radar chart data for a candidate"""
    labels = []
    values = []
    
    for metric in metrics:
        if metric in candidate:
            labels.append(metric.upper())
            # Normalize to 0-100
            val = candidate[metric]
            if isinstance(val, (int, float)):
                values.append(min(100, max(0, val * 100)))
            else:
                values.append(50)
    
    return {'labels': labels, 'values': values}

def generate_bar_chart(candidates, metric):
    """Generate bar chart data for a metric across candidates"""
    labels = [f"Candidate {c.get('rank', i+1)}" for i, c in enumerate(candidates[:10])]
    values = [c.get(metric, 0) for c in candidates[:10]]
    return {'labels': labels, 'values': values}

def create_html_table(candidates, columns):
    """Create HTML table from candidates data"""
    rows = []
    
    for cand in candidates[:50]:  # Limit to top 50
        cells = []
        for col in columns:
            val = cand.get(col, 'N/A')
            if isinstance(val, float):
                cells.append(f"<td>{val:.4f}</td>")
            else:
                cells.append(f"<td>{val}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    
    return '\n'.join(rows)

def create_html_header():
    """Generate HTML header with CSS"""
    return """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>mRNA Design Pipeline - Results Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            min-height: 100vh;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        .header p {
            opacity: 0.8;
            font-size: 1.1em;
        }
        .summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            padding: 30px;
            background: #f8f9fa;
        }
        .stat-card {
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            text-align: center;
        }
        .stat-card .value {
            font-size: 2em;
            font-weight: bold;
            color: #667eea;
        }
        .stat-card .label {
            color: #666;
            margin-top: 5px;
        }
        .section {
            padding: 30px;
        }
        .section h2 {
            color: #1a1a2e;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #667eea;
        }
        .chart-container {
            margin: 20px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 10px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        th {
            background: #1a1a2e;
            color: white;
            padding: 15px;
            text-align: left;
        }
        td {
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }
        tr:hover {
            background: #f5f5f5;
        }
        .rank-1 { background: #ffd700 !important; }
        .rank-2 { background: #c0c0c0 !important; }
        .rank-3 { background: #cd7f32 !important; }
        .badge {
            display: inline-block;
            padding: 5px 10px;
            border-radius: 20px;
            font-size: 0.8em;
            font-weight: bold;
        }
        .badge-good { background: #28a745; color: white; }
        .badge-moderate { background: #ffc107; color: #333; }
        .badge-poor { background: #dc3545; color: white; }
        .footer {
            background: #1a1a2e;
            color: white;
            padding: 20px;
            text-align: center;
            font-size: 0.9em;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 mRNA Design Pipeline</h1>
            <p>Candidate Ranking Report</p>
            <p>Generated: """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</p>
        </div>
        
        <div class="summary">
            <div class="stat-card">
                <div class="value">""" + "{total}" + """</div>
                <div class="label">Total Candidates</div>
            </div>
            <div class="stat-card">
                <div class="value">""" + "{top_score}" + """</div>
                <div class="label">Top Score</div>
            </div>
            <div class="stat-card">
                <div class="value">""" + "{avg_score}" + """</div>
                <div class="label">Average Score</div>
            </div>
            <div class="stat-card">
                <div class="value">""" + "{avg_cai}" + """</div>
                <div class="label">Average CAI</div>
            </div>
        </div>
        
        <div class="section">
            <h2>📊 Top 10 Candidates - Composite Scores</h2>
            <div class="chart-container">
                <div id="scoreChart"></div>
            </div>
        </div>
        
        <div class="section">
            <h2>📈 Score Distribution</h2>
            <div class="chart-container">
                <div id="distributionChart"></div>
            </div>
        </div>
        
        <div class="section">
            <h2>🏆 Detailed Rankings</h2>
            <table>
                <thead>
                    <tr>
                        <th>Rank</th>
                        <th>ID</th>
                        <th>Composite Score</th>
                        <th>CAI</th>
                        <th>tAI</th>
                        <th>GC%</th>
                        <th>MFE</th>
                        <th>Immunogenicity</th>
                        <th>Kozak</th>
                        <th>Sequence (first 50bp)</th>
                    </tr>
                </thead>
                <tbody>
"""

def create_html_footer():
    """Generate HTML footer with Chart.js code"""
    return """                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <p>mRNA Design Pipeline | Generated by Nextflow</p>
        </div>
    </div>
    
    <script>
        // Score chart
        var scoreData = [SCORE_DATA];
        var scoreLabels = [SCORE_LABELS];
        
        new Plotly.newPlot('scoreChart', [{
            x: scoreLabels,
            y: scoreData,
            type: 'bar',
            marker: {
                color: scoreData.map(function(v, i) {
                    return i === 0 ? '#28a745' : '#667eea';
                })
            }
        }], {
            title: 'Top 10 Candidates - Composite Scores',
            xaxis: { title: 'Candidate' },
            yaxis: { title: 'Score' }
        });
        
        // Distribution chart
        var allScores = [ALL_SCORES];
        
        new Plotly.newPlot('distributionChart', [{
            x: allScores,
            type: 'histogram',
            marker: { color: '#764ba2' }
        }], {
            title: 'Score Distribution',
            xaxis: { title: 'Composite Score' },
            yaxis: { title: 'Count' }
        });
    </script>
</body>
</html>
"""

def generate_report(ranked_data):
    """
    Generate complete HTML report
    """
    # Read ranked data from file argument
    if len(sys.argv) < 2:
        print("Usage: generate_report.py <ranked.json>", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not Path(input_file).exists():
        print(f"Error: File not found: {input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Parse JSON
    try:
        with open(input_file, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        print("Error: Invalid JSON", file=sys.stderr)
        sys.exit(1)
    
    # Extract candidates
    if isinstance(data, dict):
        candidates = data.get('all_candidates', data.get('top_10', []))
        total = data.get('total_candidates', len(candidates))
    else:
        candidates = data
        total = len(candidates)
    
    if not candidates:
        print("Error: No candidates to report", file=sys.stderr)
        sys.exit(1)
    
    # Calculate summary stats
    scores = [c.get('composite_score', 0) for c in candidates]
    top_score = max(scores) if scores else 0
    avg_score = sum(scores) / len(scores) if scores else 0
    
    cai_values = [c.get('CAI', 0) for c in candidates if 'CAI' in c]
    avg_cai = sum(cai_values) / len(cai_values) if cai_values else 0
    
    # Generate table rows
    table_rows = []
    for cand in candidates[:50]:
        rank = cand.get('rank', '-')
        rank_class = f"rank-{rank}" if rank <= 3 else ""
        
        seq_preview = cand.get('sequence', '')[:50] + '...' if len(cand.get('sequence', '')) > 50 else cand.get('sequence', '')
        
        # Immunogenicity badge
        immuno = cand.get('immunogenicity_score', 0)
        if immuno < 20:
            immuno_badge = '<span class="badge badge-good">Low</span>'
        elif immuno < 40:
            immuno_badge = '<span class="badge badge-moderate">Moderate</span>'
        else:
            immuno_badge = '<span class="badge badge-poor">High</span>'
        
        row = f"""                    <tr class="{rank_class}">
                        <td><strong>#{rank}</strong></td>
                        <td>{cand.get('id', 'N/A')}</td>
                        <td><strong>{cand.get('composite_score', 0):.2f}</strong></td>
                        <td>{cand.get('CAI', 'N/A')}</td>
                        <td>{cand.get('tAI', 'N/A')}</td>
                        <td>{cand.get('GC_content', 'N/A')}%</td>
                        <td>{cand.get('MFE', 'N/A')}</td>
                        <td>{immuno_badge} ({cand.get('immunogenicity_score', 'N/A')})</td>
                        <td>{cand.get('kozak_score', 'N/A')}</td>
                        <td><code>{seq_preview}</code></td>
                    </tr>"""
        table_rows.append(row)
    
    # Chart data
    top_10 = candidates[:10]
    score_data = ','.join([str(c.get('composite_score', 0)) for c in top_10])
    score_labels = ','.join([f"'#{c.get(\"rank\", i+1)}'" for i, c in enumerate(top_10)])
    all_scores = ','.join([str(s) for s in scores])
    
    # Build HTML
    html = create_html_header()
    html = html.replace('{total}', str(total))
    html = html.replace('{top_score}', f"{top_score:.2f}")
    html = html.replace('{avg_score}', f"{avg_score:.2f}")
    html = html.replace('{avg_cai}', f"{avg_cai:.4f}")
    
    html += '\n'.join(table_rows)
    html += create_html_footer()
    
    html = html.replace('SCORE_DATA', score_data)
    html = html.replace('SCORE_LABELS', score_labels)
    html = html.replace('ALL_SCORES', all_scores)
    
    # Output
    print(html)

def main():
    generate_report(None)

if __name__ == '__main__':
    main()
