import random

logo = """
                                  ____
                                 /\\' .\    _____
         R a n d o m            /: \___\  / .  /\\
 P h e n o t y p i c   M I M    \\' / . / /____/..\\
      G e n e r a t o r          \/___/  \\'  '\  /
                                          \\'__'\/
"""

print(logo)
print()

num_MIMs = input(' (✿◠‿◠)  Please tell me how many MIM numbers would you like to generate today:  ')

my_file = open("ALL_OMIM_MIM_NUMBERS_2022-06-10.txt", "r")
MIM_list = my_file.readlines()

print()
print('\(✿◠‿◠)/  Here you go!')
print()

for i in range(int(num_MIMs)):
    print(random.choice(MIM_list).replace('\n',''))
print()