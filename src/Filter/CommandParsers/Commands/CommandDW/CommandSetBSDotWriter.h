#ifndef SCAFFOLDER_COMMANDSETBSDORWRITER_H
#define SCAFFOLDER_COMMANDSETBSDORWRITER_H

#include <Filter/Writers/FileValidator/ValidatorNotPathWithAllLib.h>
#include <Filter/Writers/FileValidator/ValidatorFork.h>
#include <Filter/CommandParsers/State.h>
#include <Filter/CommandParsers/Commands/Command.h>
#include <Filter/Writers/DotWriter/BlockSplitDotWriterBuilder.h>


namespace filter {
    namespace commands {
        class CommandSetBSDotWriter : public Command {
        protected:

        public:
            void execute(std::string argv, State &state, ContigGraph &graph) override {
                INFO("setBlock dot writer");
                delete state.dotWriterBuilder;
                state.dotWriterBuilder = new writers::BlockSplitDotWriterBuilder();
            }
        };
    }
}

#endif //SCAFFOLDER_COMMANDSETBSDORWRITER_H
