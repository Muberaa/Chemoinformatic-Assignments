import random
import math
import matplotlib.pyplot as plt
import pandas as pd
random.seed(42) #Kodu rastgele çalıştırdığım için her seferinde rastgele aldığı parametrelerde aynı koordinatları alması için seed kullandım.
generations = 10
pop_size = 20
mutation_rate = 0.1

num_cities = 10
cities = {i: (random.uniform(0, 100), random.uniform(0, 100)) for i in range(num_cities)}

def calculate_distance(city1, city2): #iki şehir arasındaki mesafeyi öklid ile ölçüyor
    x1, y1 = cities[city1]
    x2, y2 = cities[city2]
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
def total_distance(tour): #Her bir atılan tur için toplam mesafeyi şehirlerv arası mesafeleri ve başa dönüş mesafesini eklerek hesapladım.
    return sum(calculate_distance(tour[i], tour[i + 1]) for i in range(len(tour) - 1)) + calculate_distance(tour[-1], tour[0])
def create_initial_population(pop_size): #Burada şehirleri popülasyon sırası ile özdeştirip rastgele değişen turlar (her biri bir birey) popülasyon sayısına eklenir.
    population = []
    for _ in range(pop_size):
        tour = list(cities.keys())
        random.shuffle(tour)
        population.append(tour)
    return population

def fitness(tour): #burada turun ne kadar iyi olduğunu anlamak için fitness algoritması kullandım. Fitness değeri ne kadar yüksekse turum o kadar kısa mesafeli olacaktır.
    return 1 / total_distance(tour)
def select_parents(population):
    fitness_values = [fitness(tour) for tour in population]
    total_fitness = sum(fitness_values)
    probabilities = [f / total_fitness for f in fitness_values]
    parents = random.choices(population, weights=probabilities, k=2)
    return parents
def crossover(parent1, parent2):  #her bir turu birey olarak atamıştık şimdi bu bireylerden ebeveyn seçiyoruz
    size = len(parent1)
    start, end = sorted(random.sample(range(size), 2)) #1. ebeveynnden çocuğa aktarılacak rastgele gen (şehir) seçiyoruz

    child = [None] * size
    child[start:end] = parent1[start:end] #1. ebeveynden alınan kısım ve boş kalan kısım olarak bir çocuk oluşturuyoruz

    pointer = 0
    for city in parent2:
        if city not in child:
            while child[pointer] is not None:
                pointer += 1
            child[pointer] = city #çocukta boş kalan genler için var olmayan genler 2. ebeveyn ile atanır.

    return child
def mutate(tour):  #Çocuklar için rastgele bir mutasyon oluşturur
    i, j = random.sample(range(len(tour)), 2)
    tour[i], tour[j] = tour[j], tour[i]

def genetic_algorithm(generations, pop_size, mutation_rate): #Burada artık populasyon oluşturma yeni çocukları ekleme ve mutasyon işlemlerini generasyonlar boyunca devam ettirerek bir data oluşturdum.
    population = create_initial_population(pop_size)
    best_distances = []
    generation_data = []

    for generation in range(generations):
        new_population = []

        for _ in range(pop_size // 2):
            parent1, parent2 = select_parents(population)
            child1 = crossover(parent1, parent2)
            child2 = crossover(parent2, parent1)

            if random.random() < mutation_rate:
                mutate(child1)
            if random.random() < mutation_rate:
                mutate(child2)

            new_population.extend([child1, child2])

        population = new_population

        best_tour = min(population, key=total_distance)
        best_distance = total_distance(best_tour)
        best_distances.append(best_distance)
        generation_data.append({"Generation": generation + 1, "Best Distance": best_distance, "Best Tour": best_tour})

    best_tour = min(population, key=total_distance)
    return best_tour, total_distance(best_tour), best_distances, generation_data

def plot_tour(tour, generation):
    plt.figure()
    x = [cities[city][0] for city in tour] + [cities[tour[0]][0]]
    y = [cities[city][1] for city in tour] + [cities[tour[0]][1]]
    plt.plot(x, y, marker='o', linestyle='-', color='b')
    for i, (xi, yi) in enumerate(zip(x, y)):
        plt.text(xi, yi, str(i), fontsize=12)
    plt.title(f"Generation {generation + 1}")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.show()

best_tour, best_distance, best_distances, generation_data = genetic_algorithm(generations, pop_size, mutation_rate)

results_df = pd.DataFrame(generation_data)
print("\nGeneration Results:")
print(results_df)

def plot_distances(best_distances):
    plt.figure()
    plt.plot(range(1, len(best_distances) + 1), best_distances, marker='o', linestyle='-', color='g')
    plt.title("Distance Evolution")
    plt.xlabel("Generation")
    plt.ylabel("Best Distance")
    plt.show()

plot_distances(best_distances)
