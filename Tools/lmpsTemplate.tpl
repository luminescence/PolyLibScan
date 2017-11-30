Generated datfile with some information | {{ description }}

{{ '{:> 12d}  atoms'.format(writer.count_number('particles')) }}
{% if writer.entry_is_needed('bond') -%}
	{{ '{:> 12d}  bonds'.format(writer.count_number('bonds')) }}
{%- endif %}
{% if writer.entry_is_needed('angle') -%}
	{{ '{:> 12d}  angles'.format(writer.count_number('angles')) }}
{%- endif %}
{% if writer.entry_is_needed('dihedral') -%}
	{{ '{:> 12d}  dihedrals'.format(writer.count_number('dihedrals')) }}
{%- endif %}
{{ '{:> 12d}  atom types'.format(writer.env.atom_type.values()|length) }}
{% if writer.entry_is_needed('bond') and env.bond_type.values()|length > 0 -%}
	{{'{:> 12d}  bond types'.format(env.bond_type.values()|length) }}
{%- endif %}
{% if writer.entry_is_needed('angle') and env.angle_type.values()|length > 0 -%}
    {{ '{:> 12d}  angle types'.format(env.angle_type.values()|length) }}
{%- endif %}
{% if writer.entry_is_needed('dihedral') and env.dihedral_type.values()|length > 0 -%}
    {{ '{:> 12d}  dihedral types'.format(env.dihedral_type.values()|length) }}
{%- endif %}
{% set box = env.sim_box() -%}
{{- '{:> 6.4f} {:> 6.4f} xlo xhi'.format(box[0,0], box[0,1]) }}
{{ '{:> 6.4f} {:> 6.4f} ylo yhi'.format(box[1,0], box[1,1]) }}
{{ '{:> 6.4f} {:> 6.4f} zlo zhi'.format(box[2,0], box[2,1]) }}

Masses

{% for type_ in writer.env.atom_type.values()|sort(False, False, 'Id') -%}
{{ '{:> 4d} {:> 7.2f}'.format(type_.Id,type_.mass) }}
{% endfor %}

{% for interaction in writer.env.ff.values() -%}
PairIJ Coeffs # {{ interaction.pair_type.kind }}

{% for pair in writer.pair_list(interaction) -%}
{{ writer.pair_string(pair)}}
{% endfor %}
{% endfor %}

{%- if writer.entry_is_needed('bond'): %}
Bond Coeffs

{% for bond in env.bond_type.values()|sort(False, False, 'Id') -%}
{{ env.coefficient_string(bond) }}
{% endfor -%}     			  

{%- endif -%}

{%- if writer.entry_is_needed('angle'): %}
Angle Coeffs

{% for angle in env.angle_type.values()|sort(False, False, 'Id') -%}
{{ env.coefficient_string(angle) }}
{% endfor -%}
{%- endif %}
Atoms

{% for particle in writer.chain_data('particles') -%}
{{ env.particle_str(particle)  }}
{% endfor %}
{% if writer.entry_is_needed('bond') and writer.count_number('bonds') > 0 -%}
Bonds

{% for bond in writer.chain_data('bonds') -%}
{{ '{:> 6d}{:> 4d}{:> 7d}{:> 7d}'.format(bond.Id, bond.type_.Id, bond.members[0].Id, bond.members[1].Id) }}
{% endfor -%}
{%- endif %}  
{% if writer.entry_is_needed('angle') and writer.count_number('angles') > 0 -%}
Angles

{% for angle in writer.chain_data('angles') -%}
{{ '{:> 6d}{:> 4d}{:> 7d}{:> 7d}{:> 7d}'.format(angle.Id, angle.type_.Id, 
	angle.members[0].Id, angle.members[1].Id, angle.members[2].Id) }}
{% endfor -%}
{%- endif -%}
