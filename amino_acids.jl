const AtomDict = Dict{Char,Tuple{Point,Char,Bool}}

"""
    AminoAcid(atoms, connections, offset)

atoms is a dictionary with an entry for each atom containing the location it
should be drawn, what element it is, and whether it has protons attached

connections is a list of pairs of atoms specifying all the bonds

offset describes how wide the whole amino acid will be when drawn
"""
struct AminoAcid
    atoms::AtomDict
    connections::Vector{Pair{Char,Char}}
    offset::Float64
end

function set_radius(r)
    global radius = r
    global v = 2*radius+2
    global a45 = v*sqrt(2)/2
    v30 = v*sqrt(3)/2
    h30 = v/2
    offset_n = 2*a45+2*radius+3
    offset_w = 2*v+2*radius+3
global amino_acids = Dict{Char,AminoAcid}(
'G' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false)),
                 ['N'=>'α', 'α'=>'O'],
                 offset_n),

'A' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β'],
                 offset_n),

'V' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(-h30, v+v30), 'C', true),
                           'G' => (Point(h30, v+v30), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'β'=>'G'],
                 offset_n),

'I' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(-h30, v+v30), 'C', true),
                           'G' => (Point(h30, v+v30), 'C', true),
                           'δ' => (Point(-h30, 2v+v30), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'β'=>'G', 'γ'=>'δ'],
                 offset_n),

'L' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(-h30, 2v+v30), 'C', true),
                           'D' => (Point(h30, 2v+v30), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D'],
                 offset_n),

'M' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(0, 3v), 'S', false),
                           'ϵ' => (Point(0, 4v), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'δ'=>'ϵ'],
                 offset_n),

'Y' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', false),
                           'δ' => (Point(a45, 2v+a45), 'C', true),
                           'D' => (Point(-a45, 2v+a45), 'C', true),
                           'ϵ' => (Point(a45, 3v+a45), 'C', true),
                           'E' => (Point(-a45, 3v+a45), 'C', true),
                           'ζ' => (Point(0, 3v+2a45), 'C', false),
                           'η' => (Point(0, 4v+2a45), 'O', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D', 'δ'=>'ϵ', 'D'=>'E', 'ϵ'=>'ζ', 'E'=>'ζ', 'ζ'=>'η'],
                 offset_n),

'F' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', false),
                           'δ' => (Point(a45, 2v+a45), 'C', true),
                           'D' => (Point(-a45, 2v+a45), 'C', true),
                           'ϵ' => (Point(a45, 3v+a45), 'C', true),
                           'E' => (Point(-a45, 3v+a45), 'C', true),
                           'ζ' => (Point(0, 3v+2a45), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D', 'δ'=>'ϵ', 'D'=>'E', 'ϵ'=>'ζ', 'E'=>'ζ'],
                 offset_n),

'C' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'S', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ'],
                 offset_n),

'S' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'O', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ'],
                 offset_n),

'T' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(h30, v+v30), 'O', true),
                           'G' => (Point(-h30, v+v30), 'C', true)),
               ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'β'=>'G'],
               offset_n),

'N' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(-h30, 2v+v30), 'O', false),
                           'D' => (Point(h30, 2v+v30), 'N', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D'],
                 offset_n),

'Q' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(0, 3v), 'C', false),
                           'ϵ' => (Point(-h30, 3v+v30), 'O', false),
                           'E' => (Point(h30, 3v+v30), 'N', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'δ'=>'ϵ', 'δ'=>'E'],
                 offset_n),

'D' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', false),
                           'δ' => (Point(-h30, 2v+v30), 'O', false),
                           'D' => (Point(h30, 2v+v30), 'O', false)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D'],
                 offset_n),

'E' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(0, 3v), 'C', false),
                           'ϵ' => (Point(-h30, 3v+v30), 'O', false),
                           'E' => (Point(h30, 3v+v30), 'O', false)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'δ'=>'ϵ', 'δ'=>'E'],
                 offset_n),

'R' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(0, 3v), 'C', true),
                           'ϵ' => (Point(0, 4v), 'N', true),
                           'ζ' => (Point(0, 5v), 'C', false),
                           'η' => (Point(-h30, 5v+v30), 'N', true),
                           'H' => (Point(h30, 5v+v30), 'N', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'δ'=>'ϵ', 'ϵ'=>'ζ', 'ζ'=>'η', 'ζ'=>'H'],
                 offset_n),

'K' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', true),
                           'δ' => (Point(0, 3v), 'C', true),
                           'ϵ' => (Point(0, 4v), 'C', true),
                           'ζ' => (Point(0, 5v), 'N', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'δ'=>'ϵ', 'ϵ'=>'ζ'],
                 offset_n),

'H' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(0, 0), 'C', true),
                           'O' => (Point(a45, -a45), 'C', false),
                           'β' => (Point(0, v), 'C', true),
                           'γ' => (Point(0, 2v), 'C', false),
                           'δ' => (Point(a45, 2v+v30), 'N', true),
                           'D' => (Point(-a45, 2v+v30), 'C', true),
                           'ϵ' => (Point(h30, 3v+v30), 'C', true),
                           'E' => (Point(-h30, 3v+v30), 'N', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D', 'δ'=>'ϵ', 'D'=>'E', 'ϵ'=>'E'],
                 offset_n),

'P' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', false),
                           'α' => (Point(v-a45, -a45), 'C', true),
                           'O' => (Point(2v-a45, -a45), 'C', false),
                           'β' => (Point(v-a45, v-a45), 'C', true),
                           'γ' => (Point(v-a45-h30, v+v30-a45), 'C', true),
                           'δ' => (Point(-a45, v-a45), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'δ'=>'N'],
                 offset_w),

'W' => AminoAcid(AtomDict('N' => (Point(-a45, -a45), 'N', true),
                           'α' => (Point(v-a45, -a45), 'C', true),
                           'O' => (Point(2v-a45, -a45), 'C', false),
                           'β' => (Point(-a45, v-a45), 'C', true),
                           'γ' => (Point(v-a45-h30, v+v30-a45), 'C', false),
                           'δ' => (Point(v-h30, v-a45), 'C', true),
                           'D' => (Point(v-a45-h30, 2v+v30-a45), 'C', false),
                           'ϵ' => (Point(2v-a45, v+v30-a45), 'N', true),
                           'E' => (Point(2v-a45-h30, 2v+v30-a45), 'C', false),
                           'e' => (Point(-a45, 3v+v30-a45), 'C', true),
                           'Z' => (Point(2v-a45, 3v+v30-a45), 'C', true),
                           'z' => (Point(v-a45-h30, 4v+v30-a45), 'C', true),
                           'η' => (Point(2v-a45-h30, 4v+v30-a45), 'C', true)),
                 ['N'=>'α', 'α'=>'O', 'α'=>'β', 'β'=>'γ', 'γ'=>'δ', 'γ'=>'D', 'δ'=>'ϵ', 'D'=>'E', 'D'=>'e', 'ϵ'=>'E', 'E'=>'Z', 'e'=>'z', 'Z'=>'η', 'z'=>'η'],
                 offset_w)
)
end
set_radius(1)