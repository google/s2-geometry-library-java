/*
 * Copyright 2021 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.google.common.geometry.benchmarks;

import static java.util.concurrent.TimeUnit.SECONDS;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.openjdk.jmh.results.RunResult;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.ChainedOptionsBuilder;
import org.openjdk.jmh.runner.options.CommandLineOptions;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;
import org.openjdk.jmh.runner.options.TimeValue;

/** A main() for S2 JMH benchmarks. */
public final class S2BenchmarkMain {
  private S2BenchmarkMain() {}

  public static void main(String[] argv) throws Exception {
    List<String> args = new ArrayList<>(Arrays.asList(argv));
    List<String> filteredArgs = args;

    // Command line options may override the default values for number of iterations, time per
    // iteration, etc. which are typically provided as annotations on JMH @Benchmarks. If neither an
    // annotation nor a command line option specify a value, JMH built in defaults will apply.
    CommandLineOptions cmdLineOptions = new CommandLineOptions(filteredArgs.toArray(new String[0]));
    ChainedOptionsBuilder optionsBuilder = new OptionsBuilder().parent(cmdLineOptions);

    // If the command line options didn't set forks, set it to 1 here, as the default is 5.
    if (!cmdLineOptions.getForkCount().hasValue()) {
      optionsBuilder.forks(1);
    }
    // Similarly, if a timeout wasn't set, set it to 1 minute here, as the default is 10 minutes.
    if (!cmdLineOptions.getTimeout().hasValue()) {
      optionsBuilder.timeout(new TimeValue(60, SECONDS));
    }

    Options runnerOptions = optionsBuilder.build();
    Runner runner = new Runner(runnerOptions);
    Collection<RunResult> results = runner.run();
  }
}
