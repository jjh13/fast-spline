import argparse

exp = argparse.ArgumentParser("Create an experiment")
exp.add_argument("ss_name", type=str, help="Spline Space name, check the makefile for names")
exp.add_argument("num_points", type=int, help="Number of points ")
exp.add_argument("val_width", type=float, help="Width of the window for visualization ")
exp.add_argument("--device", type=str, default="CPU")
if __name__ ==  "__main__":
    args = exp.parse_args()
    if args.device == "CPU":
        print(f"$(eval $(call SS_{args.device}_VALIDATION,{args.ss_name},{args.val_width}))")
        print(f"$(eval $(call SS_{args.device}_VEC_VALIDATION,{args.ss_name},{args.val_width}))")
    for i in range(1, args.num_points+1):

        if args.device == "GPU":
            print(f"$(eval $(call SS_{args.device}_LIN_EXPERIMENT,{args.ss_name},{i}))")

        for j in range(1, i + 1):
            if args.device == "CPU":
                print(f"$(eval $(call SS_{args.device}_VEC_EXPERIMENT,{args.ss_name},{i},{j}))")
            print(f"$(eval $(call SS_{args.device}_EXPERIMENT,{args.ss_name},{i},{j}))")
