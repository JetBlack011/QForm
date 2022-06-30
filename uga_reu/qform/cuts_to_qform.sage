import qform

cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([1, 1, 0, 1]), vector([0, 1, -1, 2])],
    [vector([1, 0, 1, 0]), vector([1, 1, 1, 1])]
]

B = qform.intersection_form(2, qform.intersection_data(2, cuts))
print(B)
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))
print()

cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([1, 2, 0, 1]), vector([0, 1, -1, 1])],
    [vector([1, 1, 1, 1]), vector([0, 1, 0, 1])]
]

B = qform.intersection_form(2, qform.intersection_data(2, cuts))
print(B)
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))


cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([1, 2, 1, 2]), vector([0, 1, 0, 1])],
    [vector([3, 2, 1, 1]), vector([1, 1, 0, 1])]
]

B = qform.intersection_form(2, qform.intersection_data(2, cuts))
print(B)
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))


cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([3, 2, 3, 2]), vector([1, 1, 1, 1])],
    [vector([1, 2, 3, 5]), vector([0, 1, 1, 2])]
]


B = qform.intersection_form(2, qform.intersection_data(2, cuts))
print(B)
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))

cuts = [
    [vector([1, 0, 0, 0]), vector([0, 1, 0, 0])],
    [vector([0, 0, 1, 0]), vector([0, 0, 0, 1])],
    [vector([3, 2, 3, 2]), vector([1, 1, 1, 1])],
    [vector([3, 2, 1, 1]), vector([1, 1, 0, 1])]
]


B = qform.intersection_form(2, qform.intersection_data(2, cuts))
print(B)
print(qform.signature(B))
print(qform.definiteness(B))
print(qform.parity(B))
