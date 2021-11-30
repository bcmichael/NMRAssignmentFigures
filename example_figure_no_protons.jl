using Luxor
radius = 3 # size of the atoms, the whole drawing is scaled to this, must be set before amino_acids.jl is included
stroke = 1 # thickness of the black lines in the drawing
include("amino_acids.jl")
include("assignments.jl")
include("drawing_utilities.jl")

atom_shapes = Dict{Char,String}('N' => "square", 'C'=>"circle", 'O'=>"triangle", 'S'=>"triangle")

seq = "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
shift_file = "example_shifts.dat"
output = "example_figure_no_protons.svg" # could be .pdf or .png etc. see Luxor documentation for all options
Drawing(540, 315, output) # set size of figure here
background("white")
sethue("black")
fontface("Helvetica")
font_size = 2*radius+2 # set the size of the text
left_start = 20 # adjust the position of the first amino acid on each row
vertical_spacing = font_size*2+5 # adjust how much vertical space is left for the labels
aa_per_row = 25
text_position = 15 # adjust how far the labels are from the drawn amino acids



fontsize(font_size)
setline(stroke)
draw_CN_legend(Point(200,5))
location = Point(left_start, vertical_spacing+3*font_size+5) # start drawing below the legend
max_height = 0
for row in Iterators.partition(assignments(seq, shift_file; M0=true), aa_per_row)
    global location, max_height
    for a in row
        location = draw_amino_acid(a..., location, radius; draw_H=false)
        max_height = max(max_height, height(a[4]))
    end
    location = Point(left_start, location.y+max_height+vertical_spacing)
    max_height = 0
end
finish()
preview()
