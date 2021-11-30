"""
    assignments(seq:String, shift_file; M0=false)

Use the sequence and a file containing a list of resonances output from sparky
to generate data structures representing all of the amino acids in the peptide
including information about which sites are assigned. Set M0 to true if the
sequence starts with a methionine at position zero.
"""
function assignments(seq::String, shift_file; M0=false)
    results = [build_assignment(aa) for aa in seq]
    shifts = split.(readlines(shift_file)[3:end])
    offset = M0 ? 1 : 0
    for shift in shifts
        shift[1] != "H2O" || continue

        code = shift[1][1]
        number = parse(Int, shift[1][2:end])+offset
        atom = uppercase(shift[2][2:end])
        element = shift[3][end]
        if atom == "" && shift[2][1] in ('N', 'H')
            atom = "N"
        end
        if (element != shift[2][1] && shift[2][1] != 'Q') || ! (atom in keys(atoms))
            continue
        end

        code == results[number].code || error()

        atom = atoms[atom]
        if element == 'H'
            results[number].protons[atom] = true
        else
            element == results[number].aa.atoms[atom][2] || println(shift)
            results[number].heavies[atom] = true
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
Translation between the atom nomenclature we use in sparky and the nomenclature
used in amino_acids.jl when describing amino acids
"""
atoms = Dict{String,Char}("A"=>'α',
                          "A2"=>'α',
                          "A3"=>'α',
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
