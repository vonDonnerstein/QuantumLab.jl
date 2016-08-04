#include <libint2.hpp>

void libint2start() {
  libint2::initialize();
}

void libint2stop() {
  libint2::finalize();
}

libint2::Shell* createShell(double* origin, int lqn, int nprim, double* exponents, double* coefficients){
  std::vector<double> alpha = std::vector<double>(exponents,exponents+nprim);
  std::vector<double> coeff = std::vector<double>(coefficients,coefficients+nprim);
  libint2::Shell* sh = new libint2::Shell(
      alpha,
      {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
      //{0, false, {0.27693435969118385, 0.26783885053249434, 0.08347366923669435}}
      {lqn, false, coeff}
      },
      {{origin[0], origin[1], origin[2]}}   // origin coordinates
      );
  return sh;
}

char lqnToLetter(int number){
  switch (number) {
    case 0: return 'S';
    case 1: return 'P'; 
    case 2: return 'D';
    case 3: return 'F';
    case 4: return 'G';
    case 5: return 'H';
    default: throw("QuantumNumber not yet implemented");
  }
}

void printShell(libint2::Shell* sh){
  printf("%c-Shell (LibInt2)   @   (%f, %f, %f)\n",lqnToLetter(sh->contr[0].l),sh->O[0],sh->O[1],sh->O[2]);

  printf("  Exponents:    ");
  for (double exp : sh->alpha) {
    printf("  %f",exp);
  }
  printf("\n");

  printf("  Coefficients: ");
  for (double exp : sh->contr[0].coeff) {
    printf("  %f",exp);
  }
  printf("\n");
}

libint2::Engine* createEngineCoulomb(int nprims, int maxlqn){
  libint2::Engine* engine= new libint2::Engine(libint2::Operator::coulomb, // will compute coulomb ints
                                               nprims,                     // max # of primitives in shells this engine will accept
                                               maxlqn,                     // max angular momentum of shells this engine will accept
                                               0
                                              );
  return engine;
}

libint2::Engine* createEngineOverlap(int nprims, int maxlqn){
  libint2::Engine* engine = new libint2::Engine(libint2::Operator::overlap,	 // will compute overlap ints
                                                                      nprims, 		 		 // max # of primitives in shells this engine will accept
                                                                      maxlqn			         // max angular momentum of shells this engine will accept
                                                                     );
  return engine;
}

libint2::Engine* createEngineKinetic(int nprims, int maxlqn){
  libint2::Engine* engine = new libint2::Engine(libint2::Operator::kinetic,	 // will compute overlap ints
                                                                      nprims, 		 		 // max # of primitives in shells this engine will accept
                                                                      maxlqn			         // max angular momentum of shells this engine will accept
                                                                     );
  return engine;
}

libint2::Engine* createEngineNuclearAttraction(int nprims, int maxlqn, libint2::Atom* atoms, int natoms){
  std::vector<libint2::Atom> atoms_vec(atoms,atoms+natoms);
  libint2::Engine* engine = new libint2::Engine(libint2::Operator::nuclear,
      						nprims,
						maxlqn
					       );
  engine->set_params(libint2::make_point_charges(atoms_vec));
  return engine;
}

void destroyEngine(libint2::Engine* engine){
  delete engine;
}
void destroyShell(libint2::Shell* sh){
  delete sh;
}

const double* compute2cInts(libint2::Engine* engine, libint2::Shell* shell1, libint2::Shell* shell2){
  engine->compute(*shell1,*shell2);
  return (engine->results())[0];
}

const double* compute4cInts(libint2::Engine* engine, libint2::Shell* mu, libint2::Shell* nu, libint2::Shell* lambda, libint2::Shell* sigma){
  engine->compute(*mu,*nu,*lambda,*sigma);
  return (engine->results())[0];
}

int main(){
  libint2::initialize();

  ////////
  printf("Create Shells.\n");
  std::vector<libint2::Shell> shells;
  shells.push_back({
	      {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.27693435969118385, 0.26783885053249434, 0.08347366923669435}}
              },
              {{0., 0., 0.}}   // origin coordinates
      });
  shells.push_back({
	      {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.27693435969118385, 0.26783885053249434, 0.08347366923669435}}
              },
              {{0., 0., 1.}}   // origin coordinates
      });
  ////////

  ////////
  printf("Create Engine.\n");
  libint2::Engine engine(libint2::Operator::overlap, // will compute overlap ints
                                        3, // max # of primitives in shells this engine will accept
                                        1      // max angular momentum of shells this engine will accept
                                       );
  ////////
  
  ////////
  printf("Compute Overlap.\n");
  for (libint2::Shell sh1 : shells) {
    for (libint2::Shell sh2 : shells) {
      engine.compute(sh1,sh2);
      printf("  %f\n",engine.results()[0][0]);
    }
  }
  ////////

  libint2::finalize();
  printf("All done.\n");
  return 0;
}
