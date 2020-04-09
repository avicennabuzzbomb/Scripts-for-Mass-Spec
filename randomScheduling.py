# random assignment of student names to each TA; demo of the class random,
# and the method pop()

import random
'''
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
'''

student_grps = ['Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5', 'Group 7', 'Group 8', 'Group 9', 'Group 10', 'Group 11', 'Group 12', 'Group 13']

TAs = ['Matt',
       'Jacob']

random.shuffle(student_grps)                    # make the list 'students' randomized
mid_index = len(student_grps) // 2              # split randomized list in half
first_half = student_grps[0:mid_index]          # separate the halves
second_half = student_grps[mid_index:]

first_TA = TAs.pop(random.randint(0, 1))    # assign each half to separate TA's
print(first_TA, sorted(first_half))
print(TAs[0], sorted(second_half))
