import qform

cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])]
]

B = qform.intersection_form(2, qform.intersection_data(2, cuts))
print(B)
print()
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))
print()
