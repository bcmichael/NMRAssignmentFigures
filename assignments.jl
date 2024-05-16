"""
    assignments(seq:String, shift_file; M0=false)

Use the sequence and a file containing a list of resonances output from sparky
to generate data structures representing all of the amino acids in the peptide
including information about which sites are assigned. Set M0 to true if the
sequence starts with a methionine at position zero.
"""
function assignments(seq::String, shift_file; M0=false)
    results = [build_assignment(aa) for aa in seq]
    shifts = readlines(shift_file)[3:end]
    offset = M0 ? 1 : 0
    for line in shifts
        shift = split(line)
        shift[1] != "H2O" || continue

        code = shift[1][1]
        number = parse(Int, shift[1][2:end])+offset
        atom = uppercase(shift[2][2:end])
        element = shift[3][end]

        # Glycine has two α protons, which might be assigned separately,
        # but will not be plotted separately here
        if code == 'G' && length(atom) == 2 && element == 'H' && atom[1] == 'A'
            atom = "A"
        end

        # The backbone N and the H attached to it do not need to
        # specify a site beyond just the element
        if atom == "" && shift[2][1] in ('N', 'H')
            atom = "N"
        end

        # The element in the atom name should match the type of nucleus
        # unless it is Q which can be used for degenerate sites
        if element != shift[2][1] && shift[2][1] != 'Q'
            println("Atom name $(shift[2]) does not match nucleus $(shift[3])")
            print_name_options(results[number].aa, element, code)
            println(line)
            continue
        end

        if ! (atom in keys(atoms))
            println("$(shift[2]) is not a valid name for an atom")
            print_name_options(results[number].aa, element, code)
            println(line)
            continue
        end

        if code != results[number].code
            println("$(shift[0]) does not match the amino acid sequence")
            println(line)
            error()
        end

        atom_name = atoms[atom]
        if ! (atom_name in keys(results[number].aa.atoms))
            println("$(shift[2]) is not a valid name for an atom in $(code)")
            print_name_options(results[number].aa, element, code)
            println(line)
            continue
        end

        if element == 'H'
            if ! results[number].aa.atoms[atom_name][3]
                site = results[number].aa.atoms[atom_name][2]*atom
                println("In $(code) $(site) does not have protons")
                print_name_options(results[number].aa, element, code)
                continue
            end
            results[number].protons[atom_name] = true
        else
            if element != results[number].aa.atoms[atom_name][2]
                println("In $(code) the element of site $(atom) is not $(element)")
                print_name_options(results[number].aa, element, code)
                println(line)
                continue
            end
            results[number].heavies[atom_name] = true
        end
    end
    return results
end

"""
    build_assignment(code::Char)

Build data structures to represent which sites are assigned in an amino acid.
The type of amino acid is speficied by its one letter code.
"""
function build_assignment(code::Char)
    aa = amino_acids[code]
    heavies = Dict{Char,Bool}()
    protons = Dict{Char,Bool}()
    for atom in aa.atoms
        heavies[atom.first] = false
        if atom.second[3] == true
            protons[atom.first] = false
        end
    end
    return (code=code, heavies=heavies, protons=protons, aa=aa)
end

"""
    print_name_options(aa::AminoAcid, element::Char, code::Char)

Print a list of the acceptable names for sites in the specified amino acid that
match the element
"""
function print_name_options(aa::AminoAcid, element::Char, code::Char)
    out = []
    for atom in atoms
        if atom[2] in keys(aa.atoms) && (aa.atoms[atom[2]][2] == element || (element == 'H' && aa.atoms[atom[2]][3]))
            push!(out, atom[1])
        end
    end
    options = element*join(out," $(element)")
    println("In $(code) valid options for $(element) are $(options)")
    # return element*join(out," $(element)")
end

"""
Translation between the atom nomenclature we use in sparky and the nomenclature
used in amino_acids.jl when describing amino acids
"""
atoms = Dict{String,Char}("A"=>'α',
                          "B"=>'β',
                          "N"=>'N',
                          "O"=>'O',
                          "G"=>'γ',
                          "G1"=>'γ',
                          "G2"=>'G',
                          "D"=>'δ',
                          "D1"=>'δ',
                          "D2"=>'D',
                          "E"=>'ϵ',
                          "E1"=>'ϵ',
                          "E2"=>'E',
                          "E3"=>'e',
                          "Z"=>'ζ',
                          "Z2"=>'Z',
                          "Z3"=>'z',
                          "H1"=>'η',
                          "H2"=>'H')
