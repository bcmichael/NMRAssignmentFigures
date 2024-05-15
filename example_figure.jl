using Luxor

include("amino_acids.jl")
include("assignments.jl")
include("drawing_utilities.jl")

atom_shapes = Dict{Char,String}('N' => "square", 'C'=>"circle", 'O'=>"triangle", 'S'=>"triangle")

seq = "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
shift_file = "example_shifts.dat"
output = "example_figure.svg" # could be .pdf or .png etc. see Luxor documentation for all options
Drawing(540, 330, output) # set size of figure here
background("white")
sethue("black")
fontface("Helvetica")
radius = 3 # size of the atoms, the whole drawing is scaled to this
stroke = 1 # thickness of the black lines in the drawing
font_size = 2*radius+2 # set the size of the text
left_start = 20 # adjust the position of the first amino acid on each row
vertical_spacing = font_size*2+5 # adjust how much vertical space is left for the labels
aa_per_row = 25
text_position = 15 # adjust how far the labels are from the drawn amino acids



fontsize(font_size)
setline(stroke)
set_radius(radius)
draw_HCN_legend(Point(200,5))
location = Point(left_start, vertical_spacing+5*font_size+5) # start drawing below the legend
max_height = 0
for row in Iterators.partition(assignments(seq, shift_file; M0=true), aa_per_row)
    global location, max_height
    for a in row
        location = draw_amino_acid(a..., location, radius)
        max_height = max(max_height, height(a[4]))
    end
    location = Point(left_start, location.y+max_height+vertical_spacing)
    max_height = 0
end
finish()
preview()
