# random assignment of student names to each TA; demo of the class random,
# and the method pop()

import random

students = ['Abdi, Lovina',
            'Busch, Samantha',
            'Cahya, Kevin',
            'Donoghue, Mary',
            'Girardi, Nicholas',
            'Gleed, Madeline',
            'Groblewski, Emma',
            'Hanson, Hallie',
            'Huang, Hui-Ching',
            'Li, Siyu',
            'Liu, Claudia',
            'Liu, Elaine',
            'Ma, Yueli',
            'Maier, Cassandra',
            'Marick, Robert',
            'Murphy, Christopher',
            'Newman, Brittany',
            'Nguyen, Thanh',
            'Schiffman, Allison',
            'Selvaraj, Monica',
            'Simmons, Zackary',
            'Sondhi, Kunal',
            'Tang, Jiayin',
            'Walls, Elijah',
            'Wermerling, Zachary']

TAs = ['Matt',
       'Jacob']

random.shuffle(students)                    # make the list 'students' randomized
mid_index = len(students) // 2              # split randomized list in half
first_half = students[0:mid_index]          # separate the halves
second_half = students[mid_index:]

first_TA = TAs.pop(random.randint(0, 1))    # assign each half to separate TA's
print(first_TA, sorted(first_half))
print(TAs[0], sorted(second_half))
