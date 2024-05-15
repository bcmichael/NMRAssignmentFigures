"""
    height(aa::AminoAcid)

Calulate how far down the amino acid comes from the origin.
"""
function height(aa::AminoAcid)
    result = 0
    for point in aa.atoms
        result = max(result, point.second[1].y)
    end
    return result+radius
end

"""
    connect_poly(p1, p2, points)

Find the point where the connecting line between p1 and p2 intersects the
polygon defined by the points.
"""
function connect_poly(p1, p2, points)
    best = [(distance(p2,p),p) for p in points[1:2]]
    if best[1][1] > best[2][1]
        reverse!(best)
    end

    for p in points[3:end]
        d = distance(p2,p)
        if d < best[2][1]
            if d < best[1][1]
                reverse!(best)
                best[1] = (d, p)
            else
                best[2] = (d, p)
            end
        end
    end
    intersectionlines(p1, p2, best[1][2], best[2][2])[2]
end

"""
    connect_point(p1, p2, element)

Find the point where the connecting line between p1 and p2 intersects the
shape for the element centered at p1.
"""
function connect_point(p1, p2, element)
    if atom_shapes[element] == "square"
        connect_poly(p1, p2, ngonside(p1, 1.75radius, 4; vertices=true))
    elseif atom_shapes[element] == "circle"
        intersectionlinecircle(p1, p2, p1, radius)[2]
    elseif atom_shapes[element] == "triangle"
        connect_poly(p1, p2, ngonside(p1, 1.75radius, 3, pi/2; vertices=true))
    else
        error()
    end
end

function draw_connecting_line(atom1, atom2)
    p1 = atom1[1]
    p2 = atom2[1]
    l1 = connect_point(p1, p2, atom1[2])
    l2 = connect_point(p2, p1, atom2[2])
    line(l1,l2, :stroke)
end

function atom_outline(atom)
    center = atom[1]
    element = atom[2]
    if atom_shapes[element] == "square"
        ngonside(center, 1.75radius, 4, 0, :stroke)
    elseif atom_shapes[element] == "circle"
        circle(center, radius, :stroke)
    elseif atom_shapes[element] == "triangle"
        ngonside(center, 1.75radius, 3, pi/2, :stroke)
    else
        error()
    end
end

function fill_atom(center, element)
    if atom_shapes[element] == "square"
        ngonside(center, 1.75radius, 4, 0, :fill)
    elseif atom_shapes[element] == "circle"
        circle(center, radius-stroke/2, :fill)
    elseif atom_shapes[element] == "triangle"
        ngonside(center, 1.75radius-1, 3, pi/2, :fill)
    else
        error()
    end
end

function draw_amino_acid(letter::Char, heavies, protons, aa::AminoAcid, position::Point, radius::Number; draw_H=true)
    origin(position)
    sethue("black")
    for n in aa.connections
        atom1 = aa.atoms[n.first]
        atom2 = aa.atoms[n.second]
        draw_connecting_line(atom1, atom2)
    end
    for n in aa.atoms
        center = n.second[1]
        element = n.second[2]
        has_protons = n.second[3] && draw_H
        atom_outline(n.second)
        assigned = false
        if element in ('C', 'N')
            color = element == 'C' ? "red" : "blue"
            if has_protons == true
                if protons[n.first] == true
                    sethue("green")
                    fill_atom(center, element)
                    assigned = true
                end
            end
            if heavies[n.first] == true
                sethue(color)
                if has_protons == true
                    circle(center, radius/2, :fill)
                else
                    fill_atom(center, element)
                    assigned = true
                end
            end
            sethue("black")
        else
            sethue("grey")
            fill_atom(center, element)
            sethue("black")
        end
    end
    sethue("black")
    text(string(letter), Point(0,-text_position), halign=:center, valign=:middle)
    CO = aa.atoms['O'][1]
    line(Point(CO.x+radius, -a45), Point(aa.offset-a45-radius,-a45))
    return Point(position.x+aa.offset, position.y)
end

function draw_HCN_legend(location)
    origin(location)
    sethue("black")
    for (x, element, color) in ((v, 'C', "red"), (2.5v, 'N', "blue"))
        text(string(element), Point(x,0), halign=:center, valign=:middle)
        atom_outline((Point(x,v), element))
        atom_outline((Point(x,2v), element))
        atom_outline((Point(x,3v), element))
        atom_outline((Point(x,4v), element))
        sethue("green")
        fill_atom(Point(x,4v), element)
        sethue(color)
        fill_atom(Point(x,2v), element)
        circle(Point(x,3v), radius/2, :fill)
        circle(Point(x,4v), radius/2, :fill)
        sethue("black")
    end
    text("O/S", Point(4v,0), halign=:center, valign=:middle)
    atom_outline((Point(4v,v), 'O'))
    sethue("grey")
    fill_atom(Point(4v,v), 'O')
    sethue("black")
    text("Unassigned", Point(0,v), halign=:right, valign=:middle)
    text("Heavy Atom Assigned (No Protons)", Point(0,2v), halign=:right, valign=:middle)
    text("Only Heavy Atom Assigned", Point(0,3v), halign=:right, valign=:middle)
    text("Heavy Atom And Protons Assigned", Point(0,4v), halign=:right, valign=:middle)
end

function draw_CN_legend(location)
    origin(location)
    sethue("black")
    for (x, element, color) in ((v, 'C', "red"), (2.5v, 'N', "blue"))
        text(string(element), Point(x,0), halign=:center, valign=:middle)
        atom_outline((Point(x,v), element))
        atom_outline((Point(x,2v), element))
        sethue(color)
        fill_atom(Point(x,2v), element)
        sethue("black")
    end
    text("O/S", Point(4v,0), halign=:center, valign=:middle)
    atom_outline((Point(4v,v), 'O'))
    sethue("grey")
    fill_atom(Point(4v,v), 'O')
    sethue("black")
    text("Unassigned", Point(0,v), halign=:right, valign=:middle)
    text("Assigned", Point(0,2v), halign=:right, valign=:middle)
end
