from riversim import *
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import time
import numpy as np

def plot(model, figsize = [40, 40], plot_limits = []):

    ax = plt.gca()

    ax.set_aspect('equal')

    plt.rcParams['figure.figsize'] = figsize
    
    ax.cla()

    xmin = 0
    xmax = 0
    ymin = 0
    ymax = 0

    for river_pr in model.rivers:
        x = ([v.x for v in river_pr.data().vertices])
        y = ([v.y for v in river_pr.data().vertices])
        ax.plot(x, y)
        if min(x) < xmin:
            xmin = min(x)
        if max(x) > xmax:
            xmax = max(x)
        if min(y) < ymin:
            ymin = min(y)
        if max(y) > ymax:
            ymax = max(y)

    
    for boundary_pr in model.region:
        x = ([v.x for v in boundary_pr.data().vertices])
        y = ([v.y for v in boundary_pr.data().vertices])
        x.append(boundary_pr.data().vertices[0].x)
        y.append(boundary_pr.data().vertices[0].y)
        if min(x) < xmin:
            xmin = min(x)
        if max(x) > xmax:
            xmax = max(x)
        if min(y) < ymin:
            ymin = min(y)
        if max(y) > ymax:
            ymax = max(y)
        ax.plot(x, y)

    
    if len(plot_limits) == 4:
        ax.set_xlim(plot_limits[0], plot_limits[1])
        ax.set_ylim(plot_limits[2], plot_limits[3])
    else:
        ax.set_xlim(xmin - 0.1, xmax + 0.1)
        ax.set_ylim(ymin - 0.1, ymax + 0.1)


    #plt.xlim(xmin - 0.1, xmax + 0.1)
    #plt.ylim(ymin - 0.1, ymax + 0.1)
    #plt.show()
    #display(fig)
    #clear_output(wait = True)
    #plt.pause(0.001)

def growRiver(m, plot_period = 50, debug = False):

    solver = Solver(m.solver_params)
    triangle = Triangle(m.mesh_params)
    mesh = TethexMesh()

    #fig = plt.figure()
    #ax = fig.add_subplot(1, 1, 1)

    dynamic_river_ids = m.rivers.tipBranchesIds()

    out = {\
        "bound_gen_t": np.zeros(m.number_of_steps), \
        "mesh_gen_t": np.zeros(m.number_of_steps), \
        "solver_t": np.zeros(m.number_of_steps), \
        "integ_t": np.zeros(m.number_of_steps), \
        "growth_t": np.zeros(m.number_of_steps), \
        "selfintersect_t": np.zeros(m.number_of_steps), \
        "plot_t": np.zeros(m.number_of_steps), \
        "total_t": np.zeros(m.number_of_steps), \
        "all_steps_time" : time.time(), \
        "status": "init"}

    for i in range(m.number_of_steps):
        out["total_t"][i] = time.time()
        if i % 10 == 0:
            print(i)
        # boundary generation: Combines boundary and river 
        # geometry into one(or several) closed boundary lines

        out["bound_gen_t"][i] = time.time()
        m.boundary = BoundaryGenerator(\
            m.sources, \
            m.region, \
            m.rivers, \
            m.region_params)
        
        num_of_intersects = NumOfBoundaryIntersection(m.boundary, m.rivers.tipBoundary())
        if num_of_intersects > 0:
            out["status"] = "self intersection detected: " + str(num_of_intersects)
            break
        out["bound_gen_t"][i] -= time.time()

        if len(m.boundary.vertices) == 0:
            out["status"] = "Boundary is empty."
            return out
    
        # mesh will be refined aroud growing tip points
        out["mesh_gen_t"][i] = time.time()
        try:
            tip_points = t_PointList()
            for id in dynamic_river_ids:
                tip_points.append(m.rivers[id].tipPoint())
            triangle.mesh_params.tip_points = tip_points
            mesh = triangle.generate(m.boundary, m.region.holes)
        except:
            out["status"] = "error in mesh generation."
            return out
        out["mesh_gen_t"][i] -= time.time()
    
        # reset solver values
        # solver.clear() lets try without it
        
        out["solver_t"][i] = time.time()
        try:
            solver.openMesh(mesh)
            for j in range(m.solver_params.adaptive_refinment_steps + 1):
                if j > 0:
                    solver.refineGrid()
                solver.setupSystem()
                solver.assembleSystem(m.boundary_conditions)
                solver.solve()
        except:
            out["status"] = "Error with solver. Most probably due to rivers intersection."
            return out
        out["solver_t"][i] -= time.time()
    
        #series parameters evaluation
        out["integ_t"][i] = time.time()
        try:
            id_series_params = t_ids_series_params()
            max_a1 = 0
            for id in dynamic_river_ids:
                tip_point = m.rivers[id].tipPoint()
                tip_angle = m.rivers[id].tipAngle()
                id_series_params[id] = solver.integrate_new(m.integr_params, tip_point, tip_angle)
                if id_series_params[id][0] > max_a1:
                    max_a1 = id_series_params[id][0]
        except:
            out["status"] = "Error in integration."
            return out
        out["integ_t"][i] -= time.time()

        if i == 100:
            print("")

        out["growth_t"][i] = time.time()
        try:
            for id_series_param in id_series_params:
                id = id_series_param.key()
                series_param = id_series_param.data()
                if m.qGrowth(series_param):
                    l = m.rivers[id].lenght()
                    if m.qBifurcate(series_param, l):
                        tip_point = m.rivers[id].tipPoint()
                        tip_angle = m.rivers[id].tipAngle()
                        br_left = Branch(tip_point, tip_angle + m.bifurcation_angle)
                        br_left.addPoint(Polar(m.ds, 0), m.region_params.river_boundary_id)
                        br_right = Branch(tip_point, tip_angle - m.bifurcation_angle)
                        br_right.addPoint(Polar(m.ds, 0), m.region_params.river_boundary_id)
                        ids = m.rivers.addSubBranches(id, br_left, br_right)

                        # add new branches
                        dynamic_river_ids.append(ids.left)
                        dynamic_river_ids.append(ids.right)

                        # and remove parent from growth evaluation
                        for f in range(len(dynamic_river_ids)):
                            if dynamic_river_ids[f] == id:
                                del dynamic_river_ids[f]
                                break
                    else: 
                        m.rivers[id].addPoint(\
                            m.nextPoint(series_param, l, max_a1),\
                            m.region_params.river_boundary_id)
                else:
                    # remove river from growth evaluation
                    for f in range(len(dynamic_river_ids)):
                        if dynamic_river_ids[f] == id:
                            del dynamic_river_ids[f]
                            break
        except:
            out["status"] = "Error in growth."
            return out
        out["growth_t"][i] -= time.time()
        
        out["plot_t"][i] = time.time()
        if i % plot_period == 0:
            plot(m)
        out["plot_t"][i] -= time.time()
        
        if debug:
            save(m, "debug " + str(i) + ".json")

        # self intersection
        out["selfintersect_t"][i] = time.time()
        empt_rivers = Rivers()
        empt_source = Sources()
        only_boundary = BoundaryGenerator(empt_source, m.region, empt_rivers, m.region_params)
        distance_to_boundary = DistanceFromPointsToBoundary(only_boundary, m.rivers.tipPoints())
        if distance_to_boundary <= m.integr_params.integration_radius:
            out["status"] = "Integration radius exceeds boundary: " + str(distance_to_boundary)
            break
        out["selfintersect_t"][i] -= time.time()

        # nothing to grow
        out["total_t"][i] -= time.time()
        if not dynamic_river_ids:
            break
    
    if out["status"] == "init":
        out["status"] = "done"

    out["all_steps_time"] -= time.time()
    
    return out