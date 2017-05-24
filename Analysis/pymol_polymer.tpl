{% for run in pymol._poly_end_poses() -%}
MODEL {{loop.index}}
{% for coords in run -%}
{{'ATOM   {:>4}  {}   LIG  {:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00            {} '.format(coords['id'], coords['ele'], coords['id'], coords['x'], coords['y'], coords['z'], coords['ele']) }}
{% endfor %}
ENDMDL
{% endfor %}