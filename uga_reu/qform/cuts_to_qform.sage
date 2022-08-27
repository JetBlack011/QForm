import qform

cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([1, 0, 1, 1]), vector([-1, 1, 0, 1])],
    [vector([1, -1, -1, -2]), vector([2, -1, 1, 1])]
]

print(qform.cuts_to_intersection(2,cuts))

B = qform.intersection_form(2, qform.cuts_to_intersection(2, cuts))
print(B)
print()
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))
print()
