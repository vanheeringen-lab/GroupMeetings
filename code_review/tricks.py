"""
General Python tips & tricks
"""
import random


"""
List comprehensions
"""
# not cool
my_list = []
for i in range(10):
    my_list.append(i)

print(my_list)

# cool
my_list = [i for i in range(10)]
print(my_list)

"""
Loop over zip
"""
list_nrs = [1, 2, 3]
list_ltr = ['a', 'b', 'c']

# not cool
for i in range(len(list_nrs)):
    nr = list_nrs[i]
    letter = list_ltr[i]
    print(nr, letter)

# cool
for nr, letter in zip(list_nrs, list_ltr):
    print(nr, letter)

"""
f-strings, since 3.6
"""
x = 123

# not cool
print("x=", x)
print("x= {}".format(x))
print("x= {var}".format(var=x))

# cool
print(f"x: {x}")
# super cool
print(f"{x=}")

"""
Type hinting, since 3.7
"""


def fibonacci(n: int) -> int:
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        return fibonacci(n - 1) + fibonacci(n - 2)


"""
Walrus operator, since 3.8
"""
# not cool
rand = random.random()
if rand < 10:
    print(f"so small: {rand=}")

# cool
if (rand := random.random()) < 0.5:
    print(f"so small: {rand=}")

"""
Classes
"""
names = ['jan', 'eva', 'piet']
ages = [10, 42, 63]
cities = ['Nijmegen', 'Amsterdam', 'Giethoorn']


class Person:
    def __init__(self, name, age, city):
        self.name = name
        self.age = age
        self.city = city

    def whoami(self):
        print(f"My name is {self.name}, I am {self.age} years old, and I live in {self.city}.")


persons = [Person(name, age, city) for name, age, city in zip(names, ages, cities)]
for person in persons:
    person.whoami()

"""
tuple unpacking
"""
# not cool
a = 'b'
b = 'a'

c = a
a = b
b = c
print(a, b)

numbers = tuple(range(100))
first_nr = numbers[0]
second_nr = numbers[1]
print(first_nr, second_nr)

# cool
a, b = 'b', 'a'
a, b = b, a
print(a, b)

numbers = tuple(range(100))
first_nr, second_nr, *_ = numbers
print(first_nr, second_nr)
