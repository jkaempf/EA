#include "evalpthread.h"

map<unsigned int, vector<double> > Evaluation::idAllelesMap;            // stock à traiter
map<unsigned int, vector<double> > Evaluation::idFitnessMap;            // stock traité

void wait(unsigned int seconds)
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}

#if ((defined(__CYGWIN__) || defined(__unix__) || defined(__MINGW32__) || defined(__MINGW64__)) && !defined(_SINGLE_THREADED))

// posix_threads
pthread_mutex_t mut1 = PTHREAD_MUTEX_INITIALIZER;

void* consumer(void *ptr) {

    Problem *problem;
    problem = (Problem *) ptr;

    unsigned int id;
    vector<double> allelesVector;

    for (;;) {

        pthread_mutex_lock( &mut1 );
        if ( Evaluation::idAllelesMap.empty() ) {
            pthread_mutex_unlock( &mut1 );
            pthread_exit(NULL); // travail terminé
        }
        else { // encore du boulot

	  id            = Evaluation::idAllelesMap.begin()->first;
	  allelesVector = Evaluation::idAllelesMap.begin()->second;

	  Evaluation::idAllelesMap.erase(Evaluation::idAllelesMap.begin());
          pthread_mutex_unlock( &mut1 );

          // evaluation proprement dite
          vector<double> fitnessVector = problem->evaluate(id, allelesVector);

          //cerr << "Recu par le thread: " << pthread_self() << "\tid" << id << "\tfitness: " << fitness << endl;

          // mettre le résultat dans la map
          pthread_mutex_lock( &mut1 );
          if ( Evaluation::idFitnessMap.size() < Evaluation::idFitnessMap.max_size() ) {
	    Evaluation::idFitnessMap.insert( pair<unsigned int,vector<double> > (id, fitnessVector) );
          }
          else throw ("Cannot allocate idFitnessMap, size full.");
          pthread_mutex_unlock( &mut1 );
        }

    }
}

#elif (defined(__QT__) && !defined(_SINGLE_THREADED))

QMutex mutex;

void Consumer::run() {

    unsigned int id;
    vector<double> allelesVector;

    for (;;) {

        mutex.lock();
        if ( Evaluation::idAllelesMap.empty() ) {
            mutex.unlock();
            return; // travail terminé
        }
        else { // encore du boulot

            id            = Evaluation::idAllelesMap.begin()->first;
            allelesVector = Evaluation::idAllelesMap.begin()->second;

            Evaluation::idAllelesMap.erase(idAllelesMap.begin());
            mutex.unlock();

            // evaluation proprement dite
            //mutex.lock();
            vector<double> fitnessVector = problem->evaluate(id, allelesVector);
            //mutex.unlock();

            //cerr << "Recu par le thread: " << pthread_self() << "\tid" << id << "\tfitness: " << fitness << endl;

            // mettre le résultat dans la map
            mutex.lock();
            if ( Evaluation::idFitnessMap.size() < Evaluation::idFitnessMap.max_size() ) {
                Evaluation::idFitnessMap.insert( pair<unsigned int,vector<double> > (id, fitnessVector) );
            }
            else throw ("Cannot allocate idFitnessMap, size full.");
            mutex.unlock();
        }

    }
}

#endif
