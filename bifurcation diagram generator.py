import numpy as np
import matplotlib.pyplot as plt
import numba

@numba.jit(nopython=True, parallel=True)
def calc(func):
  # je nachdem welche Funktion ausgewählt wurde, haben wir hier
  #die Charakterisierung des Parametergitters
  #max_iterations lässt sich verändern um die Qualität der Bifurkationsdiagramme zu verbessern
  match func:
    case 1:
        min_r = 2.4
        max_r = 4
        step_r = 0.001
        max_iterations = 1000
    case 2:
        min_r = 1
        max_r = 2
        step_r = 0.0005
        max_iterations = 900
    case 3:
        min_r = 2
        max_r = 3
        step_r = 0.0005
        max_iterations = 1000
    case 4:
        min_r = 3.5
        max_r = 4
        step_r = 0.0002
        max_iterations = 700
    case 5:
        min_r = 1
        max_r = 2
        step_r = 0.0005
        max_iterations = 700
  skip_iterations = 500
  max_count = int((max_iterations - skip_iterations) * (max_r - min_r) / step_r)
  result_x = np.zeros(max_count)
  result_r = np.zeros(max_count)
  i = 0
  for r in np.arange(min_r, max_r, step_r):
    x = 0.1
    for it in range(max_iterations):
        match func:
            case 1:
                x = r * x * (1-x)
            case 2:
                x = r * min(x,(1-x))
            case 3:
                x = (1/(1-pow(2,1-r)))*(1-pow(x,r)-pow((1-x),r))
            case 4:
                x = r * x * (1-x)
            case 5:
                x = r * min(x,(1-x))
        if it > skip_iterations:
            result_x[i] = x
            result_r[i] = r
            i += 1
  result_x = result_x[result_r != 0].copy()
  result_r = result_r[result_r != 0].copy()
  return result_x, result_r


@numba.jit(nopython=True, parallel=True)
def isolines(func,max_iso):#erstellt die ersten max_iso grenzfunktionen unserer Funktion (log oder tentmap)
  # Parameters
  match func:
    case 1:
      min_r = 3.5
      max_r = 4
      step_r = 0.0002
    case 2:
      min_r = 1
      max_r = 2
      step_r = 0.0005
      max_iterations = 700
  max_counter = int( (max_r - min_r) / step_r)
  result_iso = np.zeros((max_iso,max_counter))
  result_r = np.zeros(max_counter)
  i = 0
  for r in np.arange(min_r, max_r, step_r):
    result_r[i] = r
    x = 0.5
    for it in range(max_iso):
      match func:
            case 1:
                x = r * x * (1-x)
            case 2:
                x = r * min(x,(1-x))
      result_iso[it][i] = x
    i += 1
  return result_iso, result_r




result_x, result_r = calc(1)

# Plot
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(result_r, result_x, ",", color='k')
plt.xlim([2.4, 4])
plt.xlabel("θ")
plt.ylabel("x")
plt.show()

result_x, result_r = calc(2)

# Plot
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(result_r, result_x, ",", color='k')
plt.xlim([1, 2])
plt.xlabel("θ")
plt.ylabel("x")

plt.show()


result_x, result_r = calc(3)
# Plot
plt.figure(figsize=(10, 6), dpi=200)
plt.plot(result_r, result_x, ",", color='k')
plt.xlim([2, 3])
plt.xlabel("θ")
plt.ylabel("x")
plt.show()


result_x, result_r = calc(4)

# Plot
plt.figure(figsize=(6, 6), dpi=200)
plt.plot(result_r, result_x, ",", color='gray')
max_isolines=8
result_iso, result_r = isolines(1,max_isolines)
for it in range(max_isolines):
    plt.plot(result_r, result_iso[it], "-",linewidth='0.3', color='k')
    
plt.xlim([3.5, 4])
plt.xlabel("θ")
plt.ylabel("x")
plt.show()


result_x, result_r = calc(5)

# Plot
plt.figure(figsize=(6, 6), dpi=200)
plt.plot(result_r, result_x, ",", color='gray')
max_isolines=8
result_iso, result_r = isolines(2,max_isolines)
for it in range(max_isolines):
    plt.plot(result_r, result_iso[it], "-",linewidth='0.3', color='k')
    
plt.xlim([1, 2])
plt.xlabel("θ")
plt.ylabel("x")
plt.show()


